library(plyr)
library(dplyr)
library(edgeR)
library(epiDisplay)
library(rsq)
library(MASS)
library(Biocomb)
library(caret)
library(dplyr)
library(epiDisplay)
library(pROC)
library(ggplot2)
library(DMwR)
library(ROSE)
library(gridExtra)

# DE funkcja
ks.miRNA.de = function(ttpm_pofiltrze, klasy)
{
  wyniki = data.frame(miR = colnames(ttpm_pofiltrze))

  library(dplyr)

  for (i in 1:length(colnames(ttpm_pofiltrze)))
  {
    # Filtry
    #wyniki[i,"jakieszero"] = ifelse(sum(ttpm_pofiltrze[,i] == 0) > 0, "tak", "nie")

    # Średnia i SD
    wyniki[i,"mean logtpm"] = mean(ttpm_pofiltrze[,i])
    wyniki[i,"median logtpm"] = median(ttpm_pofiltrze[,i])
    wyniki[i,"SD logtpm"] = sd(ttpm_pofiltrze[,i])

    # Cancer
    tempx = ttpm_pofiltrze[klasy == "Cancer",]
    wyniki[i,"cancer mean"] = mean(tempx[,i])
    wyniki[i,"cancer median"] = median(tempx[,i])
    wyniki[i,"cancer SD"] = sd(tempx[,i])

    # Cancer
    tempx = ttpm_pofiltrze[klasy != "Cancer",]
    wyniki[i,"control mean"] = mean(tempx[,i])
    wyniki[i,"control median"] = median(tempx[,i])
    wyniki[i,"control SD"] = sd(tempx[,i])

    # DE
    temp = t.test(ttpm_pofiltrze[,i] ~ as.factor(klasy))
    fc = (temp$estimate[1] - temp$estimate[2])
    wyniki[i,"log10FC"] = fc
    wyniki[i,"log2FC"] = wyniki[i,"log10FC"] / log10(2)

    revfc = (temp$estimate[2] - temp$estimate[1])
    wyniki[i,"reversed_log10FC"] = revfc
    wyniki[i,"reverse_log2FC"] = wyniki[i,"reversed_log10FC"] / log10(2)
    #wyniki[i,"log2FC"] = log2(fc)
    wyniki[i,"p-value"] = temp$p.value
  }
  wyniki[,"p-value Bonferroni"] = p.adjust(wyniki$`p-value`, method = "bonferroni")
  wyniki[,"-log10(p-value Bonferroni)"] = -log10(p.adjust(wyniki$`p-value`, method = "bonferroni"))
  wyniki[,"p-value BH"] = p.adjust(wyniki$`p-value`, method = "BH")
  return(wyniki)
}

ks.wykreskorelacji = function(var1, var2, labvar1, labvar2, title) {
  plot(var1, var2, pch = 19, cex=0.5, xlab=labvar1, ylab=labvar2, main=title)
  abline(0,1, col='gray')
  abline(fit <- lm(var2 ~ var1), col='black', lty = 2)
  temp = cor.test(var1, var2, method = 'pearson')
  legend("topleft", bty="n", legend=paste("r =",round(cor(var1,var2),2),", p =",round(temp$p.value, 4),"\nR² =", format(summary(fit)$adj.r.squared, digits=4)))
}

ks.miR.formula = function(wybrane_miRy) {
  as.formula(paste0("Class ~ ",paste0(as.character(wybrane_miRy), collapse = " + ")))
}

ks.wczytajmix = function(wd = getwd(), smote_over = 10000, use_smote_not_rose = F, replace_smote = F) {
  oldwd = getwd()
  setwd(wd)
  train = dplyr::select(read.csv("mixed_train.csv", stringsAsFactors = F), starts_with("hsa"), Class)

  test = dplyr::select(read.csv("mixed_test.csv", stringsAsFactors = F), starts_with("hsa"), Class)
  valid = dplyr::select(read.csv("mixed_valid.csv", stringsAsFactors = F), starts_with("hsa"), Class)
  train$Class = factor(train$Class, levels = c("Control","Cancer"))
  test$Class = factor(test$Class, levels = c("Control","Cancer"))
  valid$Class = factor(valid$Class, levels = c("Control","Cancer"))

  # Wywalamy miRy z zerowa wariancja
  temp = train %>% filter(Class == "Cancer")
  temp2 = as.numeric(which(apply(temp, 2, var) == 0))
  temp = train %>% filter(Class == "Control")
  temp3 = as.numeric(which(apply(temp, 2, var) == 0))
  temp4 = unique(c(temp2, temp3))
  if (length(temp4) > 0) {
  train = train[,-temp4]
  test = test[,-temp4]
  valid = valid[,-temp4]
  }

  if(use_smote_not_rose) {
    train_smoted = DMwR::SMOTE(Class ~ ., data = train, perc.over = smote_over,perc.under=100, k=10)
  train_smoted$Class = factor(train_smoted$Class, levels = c("Control","Cancer"))
  } else {
  rosed = ROSE(Class ~ ., data = train, N = nrow(train)*10, seed = 1)
  train_smoted = rosed[["data"]]
  train_smoted$Class = factor(train_smoted$Class, levels = c("Control","Cancer"))
  }



  if (replace_smote == T) { train = train_smoted }

  trainx = dplyr::select(train, starts_with("hsa"))
  trainx_smoted = dplyr::select(train_smoted, starts_with("hsa"))
  setwd(oldwd)
  return(list(train, test, valid, train_smoted, trainx, trainx_smoted))
}

ks.selekcjamiRNA = function(wd = getwd(), m = c(1:50),
                            max_iterations = 10, code_path = "/home/konrad/snorlax/2019_PRELUDIUM/feature_selection/",
                            register_parallel = T, clx = NULL, stamp = as.numeric(Sys.time())) {

  oldwd = getwd()
  setwd(wd)


  if(!dir.exists("temp")) { dir.create("temp") }

  run_id = stamp
  formulas = list()
  times = list()

  zz <- file(paste0("temp/",stamp,"featureselection.log"), open = "wt")
  sink(zz)
  #sink(zz, type = "message")
  pdf(paste0("temp/",stamp,"featureselection.pdf"))




  dane = ks.wczytajmix(); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  n = 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
  start_time <- Sys.time()

  formulas[["all"]] = ks.miR.formula(colnames(trainx))

  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  cat("\n\nStandard DE...")
  wyniki = ks.miRNA.de(trainx, train$Class)
  istotne = filter(wyniki, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)

  train_sig = dplyr::select(train, as.character(istotne$miR), Class)
  trainx_sig = dplyr::select(train, as.character(istotne$miR))

  cat("\n\nPreparing for SMOTE...")
  #train_sig_smoted = DMwR::SMOTE(Class ~ ., data = train_sig, perc.over = 10000,perc.under=100, k=10)
  train_sig_smoted = dplyr::select(train_smoted, as.character(istotne$miR), Class)
  train_sig_smoted$Class = factor(train_sig_smoted$Class, levels = c("Control","Cancer"))
  trainx_sig_smoted = dplyr::select(train_sig_smoted, starts_with("hsa"))

  wyniki_smoted = ks.miRNA.de(trainx_smoted, train_smoted$Class)
  istotne_smoted = filter(wyniki_smoted, `p-value BH` <= 0.05) %>% arrange(`p-value BH`)

  # Caret prep
  if (register_parallel) {
  cat("\n\nGetting cluster ready...")
  if(is.null(clx)) {
  library(doParallel)
  cl <- makePSOCKcluster(detectCores() - 2)
  registerDoParallel(cl) }
    else { registerDoParallel(clx) }
  }

  # 0. All and sig
  cat("\n\nStarting selected methods...")


  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
    cat("\n\nStarting SIG")
    start_time <- Sys.time()

  formulas[["sig"]] = ks.miR.formula(as.character(istotne$miR))


  formulas[["sigSMOTE"]] = ks.miR.formula(as.character(istotne$miR))

  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  # 1. Fold-change and sig. filter
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
  start_time <- Sys.time()

  cat("\n\nStarting FC and SIG")
  fcsig = as.character(istotne$miR[abs(istotne$log2FC)>1])
  formulas[["fcsig"]] = ks.miR.formula(fcsig)

  fcsig = as.character(istotne_smoted$miR[abs(istotne_smoted$log2FC)>1])
  formulas[["fcsigSMOTE"]] = ks.miR.formula(fcsig)

  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  # 2. CFS
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
    start_time <- Sys.time()
  cat("\n\nStarting CFS")
  cfs = select.cfs(train)
  formulas[["cfs"]] = ks.miR.formula(as.character(cfs$Biomarker))

  cfs = select.cfs(train_smoted)
  formulas[["cfsSMOTE"]] = ks.miR.formula(as.character(cfs$Biomarker))

  cfs = select.cfs(train_sig)
  formulas[["cfs_sig"]] = ks.miR.formula(as.character(cfs$Biomarker))

  cfs = select.cfs(train_sig_smoted)
  formulas[["cfsSMOTE_sig"]] = ks.miR.formula(as.character(cfs$Biomarker))
  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  # 3. classifier.loop
  n = n + 1
  if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting.."));
  start_time <- Sys.time()
  cat("\n\nStarting classifier loop")
  classloop = classifier.loop(train, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=10)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloop"]] = ks.miR.formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_smoted, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=10)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloopSMOTE"]] = ks.miR.formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_sig, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=10)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloop_sig"]] = ks.miR.formula(f_classloop)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classloop = classifier.loop(train_sig, feature.selection = "auc", method.cross="fold-crossval", classifiers=c("svm","lda","rf","nsc"), no.feat=10)
  f_classloop = rownames(classloop$no.selected)[classloop$no.selected[,1]>0]
  formulas[["classloopSMOTE_sig"]] = ks.miR.formula(f_classloop)

  end_time <- Sys.time()
  saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS"))
  saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  # 4. select.forward.Corr
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting select.forward.Corr")
  fcfs = select.forward.Corr(train, disc.method="MDL")
  formulas[["fcfs"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_smoted, disc.method="MDL")
  formulas[["fcfsSMOTE"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_sig, disc.method="MDL")
  formulas[["fcfs_sig"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fcfs = select.forward.Corr(train_sig_smoted, disc.method="MDL")
  formulas[["fcfsSMOTE_sig"]] = ks.miR.formula(fcfs)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # 5. select.forward.wrapper
  fwrap = select.forward.wrapper(train)
  formulas[["fwrap"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_smoted)
  formulas[["fwrapSMOTE"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_sig)
  formulas[["fwrap_sig"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  fwrap = select.forward.wrapper(train_sig_smoted)
  formulas[["fwrapSMOTE_sig"]] = ks.miR.formula(fwrap)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  # 6, 7, 8.
  #select.process(dattable,method="InformationGain",disc.method="MDL",
  #               threshold=0.2,threshold.consis=0.05,attrs.nominal=numeric(),
  #               max.no.features=10)
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting MDL")
  formulas[["AUC_MDL"]] = ks.miR.formula(colnames(train)[select.process(train, method="auc", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDL"]] = ks.miR.formula(colnames(train)[select.process(train, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDL"]] = ks.miR.formula(colnames(train)[select.process(train, method="CorrSF", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();

  formulas[["AUC_MDLSMOTE"]] = ks.miR.formula(colnames(train_smoted)[select.process(train_smoted, method="auc", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();


  formulas[["SU_MDLSMOTE"]] = ks.miR.formula(colnames(train_smoted)[select.process(train_smoted, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDLSMOTE"]] = ks.miR.formula(colnames(train_smoted)[select.process(train_smoted, method="CorrSF", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["AUC_MDL_sig"]] = ks.miR.formula(colnames(train_sig)[select.process(train_sig, method="auc", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDL_sig"]] = ks.miR.formula(colnames(train_sig)[select.process(train_sig, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = 10)])

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDL_sig"]] = ks.miR.formula(colnames(train_sig)[select.process(train_sig, method="CorrSF", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["AUC_MDLSMOTE_sig"]] = ks.miR.formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="auc", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["SU_MDLSMOTE_sig"]] = ks.miR.formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="symmetrical.uncertainty", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  formulas[["CorrSF_MDLSMOTE_sig"]] = ks.miR.formula(colnames(train_sig_smoted)[select.process(train_sig_smoted, method="CorrSF", disc.method = "MDL", max.no.features = 10)])
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  # 9. bounceR - genetic
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting bounceR..")
  library(bounceR)
  mrmr <- featureSelection(data = train,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = 10,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full"]] = mrmr@opt_formula
  formulas[["bounceR-stability"]] = ks.miR.formula(as.character(mrmr@stability[1:10,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  mrmr <- featureSelection(data = train_smoted,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = 10,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SMOTE"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SMOTE"]] = ks.miR.formula(as.character(mrmr@stability[1:10,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();

  mrmr <- featureSelection(data = train_sig,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = 10,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SIG"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SMOTE"]] = ks.miR.formula(as.character(mrmr@stability[1:10,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();

  mrmr <- featureSelection(data = train_sig_smoted,
                           target = "Class",
                           max_time = "15 mins",
                           selection = selectionControl(n_rounds = NULL,
                                                        n_mods = NULL,
                                                        p = 10,
                                                        penalty = 0.5,
                                                        reward = 0.2),
                           bootstrap = "regular",
                           early_stopping = "aic",
                           boosting = boostingControl(mstop = 100, nu = 0.1),
                           cores = parallel::detectCores()-1)
  formulas[["bounceR-full_SIGSMOTE"]] = mrmr@opt_formula
  formulas[["bounceR-stability_SIGSMOTE"]] = ks.miR.formula(as.character(mrmr@stability[1:10,] %>% pull('feature')))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }


  # 10. RFE RandomForest
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting RFE RF")
  ctrl <- rfeControl(functions =rfFuncs,
                     method = "cv", number = 10,
                     saveDetails = TRUE,
                     allowParallel = TRUE,
                     returnResamp = "all",
                     verbose = T)

  rfProfile <- rfe(trainx, train$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFE"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rfProfile <- rfe(trainx_smoted, train_smoted$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFESMOTE"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rfProfile <- rfe(trainx_sig, train_sig$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFE_sig"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rfProfile <- rfe(trainx_sig_smoted, train_sig_smoted$Class,
                   sizes = 3:11,
                   rfeControl = ctrl)
  plot(rfProfile, type=c("g", "o"))
  formulas[["RandomForestRFESMOTE_sig"]] = ks.miR.formula(predictors(rfProfile))
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }



  # 11. Genetic
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting Genetic")
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx, y = train$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRF"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_smoted, y = train_smoted$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRFSMOTE"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_sig, y = train_sig$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRF_sig"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  ga_ctrl <- gafsControl(functions = rfGA, method = "repeatedcv", number=10, repeats=5, allowParallel=T)
  rf_ga <- gafs(x = trainx_sig_smoted, y = train_sig_smoted$Class,
                iters = max_iterations,
                gafsControl = ga_ctrl)
  plot(rf_ga) + theme_bw()
  print(rf_ga)
  formulas[["GeneticAlgorithmRFSMOTE_sig"]] = ks.miR.formula(rf_ga$ga$final)
  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
 save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }


  # 12. SimulatedAnealing
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting SimulatedAnealing")
  sa_ctrl <- safsControl(functions = rfSA,
                         method = "repeatedcv",
                         number=5, repeats=10, allowParallel=T,
                         improve = 50)

  rf_sa <- safs(x = trainx, y = train$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRF"]] = ks.miR.formula(rf_sa$sa$final)



  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rf_sa <- safs(x = trainx_smoted, y = train_smoted$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRFSMOTE"]] = ks.miR.formula(rf_sa$sa$final)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rf_sa <- safs(x = trainx_sig, y = train_sig$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRF_sig"]] = ks.miR.formula(rf_sa$sa$final)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  rf_sa <- safs(x = trainx_sig_smoted, y = train_sig_smoted$Class,
                iters = max_iterations,
                safsControl = sa_ctrl)
  print(rf_sa)
  plot(rf_sa) + theme_bw()
  formulas[["SimulatedAnnealingRFSMOTE_sig"]] = ks.miR.formula(rf_sa$sa$final)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }



  # 14. Boruta (https://www.jstatsoft.org/article/view/v036i11)
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  library(Boruta)
  cat("\n\nStarting Boruta")
  bor = Boruta(trainx, train$Class)
  formulas[["Boruta"]] = ks.miR.formula(names(bor$finalDecision)[as.character(bor$finalDecision) == "Confirmed"])

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  bor = Boruta(trainx_smoted, train_smoted$Class)
  formulas[["BorutaSMOTE"]] = ks.miR.formula(names(bor$finalDecision)[as.character(bor$finalDecision) == "Confirmed"])

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }


  # 15. spFSR - feature selection and ranking by simultaneous perturbation stochastic approximation
  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  cat("\n\nStarting spFSR")
  library(spFSR)
  knnWrapper    <- makeLearner("classif.knn", k = 5)
  classifTask   <- makeClassifTask(data = train, target = "Class")
  perf.measure  <- acc
  spsaMod <- spFeatureSelection(
    task = classifTask,
    wrapper = knnWrapper,
    measure = perf.measure ,
    num.features.selected = 10,
    iters.max = max_iterations,
    num.cores = detectCores() - 1)
  formulas[["spFSR"]] = ks.miR.formula(spsaMod$features)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  classifTask   <- makeClassifTask(data = train_smoted, target = "Class")
  perf.measure  <- acc
  spsaMod <- spFeatureSelection(
    task = classifTask,
    wrapper = knnWrapper,
    measure = perf.measure ,
    num.features.selected = 10,
    iters.max = max_iterations,
    num.cores = detectCores() - 2)
  formulas[["spFSRSMOTE"]] = ks.miR.formula(spsaMod$features)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # varSelRF
  cat("\n\nStarting varSelRF")
  library(varSelRF)
  var.sel <- varSelRF(trainx, train$Class, ntree = 500, ntreeIterat = max_iterations, vars.drop.frac = 0.05, whole.range = T, keep.forest = T)
  formulas[["varSelRF"]] = ks.miR.formula(var.sel$selected.vars)

  var.sel <- varSelRF(trainx_smoted, train_smoted$Class, ntree = 500, ntreeIterat = max_iterations, vars.drop.frac = 0.05, whole.range = T, keep.forest = T)
  formulas[["varSelRFSMOTE"]] = ks.miR.formula(var.sel$selected.vars)

  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
  # 13. WxNet (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6642261/)
  cat("\nStarting WxNet")
  library(plyr)
  library(dplyr)
  library(reticulate)
  library(tidyverse)
  library(data.table)
  library(DMwR)



  # Set the path to the Python executale file
  #use_python("/anaconda3/bin/python", required = T)


  conda_list()
  use_condaenv("tensorflow", required = T)
  py_config()

  dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  #train_org = train
  #train = SMOTE(Class ~ ., data = train_org, perc.over = 10000,perc.under=100)

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX

  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )





  # Wx with SMOTE
  dane = ks.wczytajmix(replace_smote = T); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta.py", local = T)
  #
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["WxSMOTE"]] = ks.miR.formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["WxSMOTE"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )

  dane = ks.wczytajmix(replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta_z.py", local = T)
  #
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["Wx_Zscore"]] = ks.miR.formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta_z.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx_Zscore"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )

  dane = ks.wczytajmix(replace_smote = T); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  # Przygotowanie
  trainx = dplyr::select(train, starts_with("hsa"))
  testx = dplyr::select(test, starts_with("hsa"))
  trainx = t(trainx)
  testx = t(testx)
  ids = paste0("train",1:ncol(trainx))
  ids2 = paste0("test",1:ncol(testx))
  colnames(trainx) = gsub("-", "", ids)
  colnames(testx) = gsub("-", "", ids2)

  traindata = cbind(data.frame(fnames = rownames(trainx), trainx))
  fwrite(traindata, paste0(code_path, "wx/DearWXpub/src/train-data.csv"), row.names = F, quote = F)
  testdata = cbind(data.frame(fnames = rownames(testx), testx))
  fwrite(testdata, paste0(code_path, "wx/DearWXpub/src/test-data.csv"), row.names = F, quote = F)

  trainanno = data.frame(id = colnames(trainx), label = train$Class)
  fwrite(trainanno, paste0(code_path, "wx/DearWXpub/src/train-anno.csv"), row.names = F, quote = F)
  testanno = data.frame(id = colnames(testx), label = test$Class)
  fwrite(testanno, paste0(code_path, "wx/DearWXpub/src/test-anno.csv"), row.names = F, quote = F)


  # Wywolanie WX
  # setwd(paste0(code_path,"wx/DearWXpub/src/"))
  # py_run_file("wx_konsta_z.py", local = T)
  #
  # np <- import("numpy")
  # w = np$load("wyniki.npy",allow_pickle=T)
  # print(w)
  # setwd(wd)
  # formulas[["Wx_ZscoreSMOTE"]] = ks.miR.formula(w)
  out <- tryCatch(
    {
      setwd(paste0(code_path,"wx/DearWXpub/src/"))
      py_run_file("wx_konsta_z.py", local = T)

      np <- import("numpy")
      w = np$load("wyniki.npy",allow_pickle=T)
      formulas[["Wx_ZscoreSMOTE"]] = ks.miR.formula(w)
    },
    error=function(cond) {
      message("ERROR:")
      message(cond)
      # Choose a return value in case of error
    },
    warning=function(cond) {
      message("WARNING:")
      message(cond)
    },
    finally={
      setwd(wd)
    }
  )


  end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
  save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }

  ######################################################### MK  ######################################################
  ## MK feature selection using multiple iterations of RFE (on radnom forrest) with 'votin' for mirs on each

  iteratedRFE <- function(trainSet, testSet, initFeatures = colnames(trainSet), classLab, checkNFeatures = 25, votingIterations = 100000, useCV = F, nfolds = 10, initRandomState = 42 ) {

    set.seed(initRandomState)

    #prepare output data structures
    resAcc <- c(rep(0, checkNFeatures))
    resVotes <- data.frame(matrix(0, nrow = length(initFeatures), ncol = checkNFeatures), row.names = initFeatures)
    for(i in 1:checkNFeatures) colnames(resVotes)[i] <- toString(i)
    resTop <- list()


    for (i in 1:votingIterations) {

      if(useCV == F) {
        params <- rfeControl(functions = rfFuncs, saveDetails = T)
        iter <- rfeIter(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), testX = testSet[, initFeatures], testY = as.factor(testSet[, classLab]), sizes = 1:checkNFeatures,
                        metric = "Accuracy", rfeControl = params)

        for(j in 1:checkNFeatures) {
          tmp <- iter$pred[iter$pred$Variables == j, ]

          acc <- length(which(tmp$pred == tmp$obs)) / nrow(tmp) #calculate and add accuracy
          resAcc[j] <- resAcc[j] + acc

          selected <- iter$finalVariables[[j+1]]$var
          numb <- iter$finalVariables[[j+1 ]]$Variables[1]

          resVotes[selected, numb] <- resVotes[selected, numb] + 1
        }


      }
      else {

        seeds <- vector(mode = "list", length = nfolds + 1) # add random seeds for cross validation
        for(i in 1:nfolds) seeds[[i]] <- sample.int(1000000000, checkNFeatures + 1)
        seeds[nfolds + 1] <- sample.int(1000000000, 1)

        params <- rfeControl(functions = rfFuncs, number = nfolds, saveDetails = T)
        iter <- rfe(x = trainSet[, initFeatures], y = as.factor(trainSet[, classLab]), sizes = 1:checkNFeatures, rfeControl = params)

        for(j in 1:checkNFeatures) {
          tmp <- iter$variables[iter$variables$Variables == j, ]
          for(k in tmp$var) resVotes[k, j] <- resVotes[k, j] + 1 # increase a voting score for each fold

          resAcc[j] <- resAcc[j] + iter$results[iter$results$Variables == j, "Accuracy"]

        }
      }
    }

    resAcc <- resAcc / votingIterations #make average accuracy

    for(i in 1:ncol(resVotes)) resTop[[i]] <- rownames(resVotes[order(-resVotes[, i])[1:i], ])

    returning <- list(data.frame(resAcc), resVotes, resTop)
    names(returning) <- c("accuracyPerNFeatures", "votesPerN", "topFeaturesPerN")
    return(returning)

  }



  n = n + 1; if (n %in% m) { cat(paste0("\n\nMatched method ", n, " with those requested.. Starting..")); start_time <- Sys.time();
    cat("\n\nStarting multiply iterated RFE")

    selectedMirsCV <- iteratedRFE(trainSet = train, classLab = 'Class', checkNFeatures = ?How many feat?)$topFeaturesPerN[[?How many feat?]]
    selectedMirsTest <- iteratedRFE(trainSet = train, testSet = test, classLab = 'Class', checkNFeatures = ?How many feat?)$topFeaturesPerN[[?How many feat?]]

    formulas[["iteratedRFECV"]] = ks.miR.formula(selectedMirsCV$topFeaturesPerN[[?How many feat?]])
    formulas[["iteratedRFETest"]] = ks.miR.formula(selectedMirsTest$topFeaturesPerN[[?How many feat?]])

    end_time <- Sys.time(); saveRDS(end_time - start_time, paste0("temp/",n,"-",run_id,".RDS")); saveRDS(formulas, paste0("temp/formulas",run_id,"-",n,".RDS"))
    save(list = ls(), file = paste0("temp/all",run_id,".rdata")); print(formulas)
  }
  ######################################################### end of MK ######################################################
  # End
  cat("\n\nEnding...")
  stopCluster(cl)
  saveRDS(formulas, paste0("temp/featureselection_formulas",stamp,".RDS"))
  saveRDS(formulas, paste0("featureselection_formulas.RDS"))
  nazwa = paste0("temp/featureselection_data",stamp,".rda")
  save(list = ls(), file = nazwa)
  dev.off()
  print(formulas)
  setwd(oldwd)
  return(formulas)
}

ks.benchmark = function(wd = getwd(), search_iters = 200, cores = detectCores()-1, input_formulas = readRDS("featureselection_formulas_final.RDS"),
                        output_file = "benchmark.csv", mxnet = F, algorithms = c("rf","C5.0", "nnet","pcaNNet"), holdout = T, stamp = as.character(as.numeric(Sys.time()))) {
  library(plyr)
  library(dplyr)
  library(edgeR)
  library(epiDisplay)
  library(rsq)
  library(MASS)
  library(Biocomb)
  library(caret)
  library(dplyr)
  library(epiDisplay)
  library(pROC)
  library(ggplot2)
  library(DMwR)
  library(stringr)
  library(psych)
  library(C50)
  library(randomForest)
  library(nnet)
  library(reticulate)
  library(stargazer)

  if(!dir.exists("temp")) { dir.create("temp") }

  use_condaenv("tensorflow")
  zz <- file(paste0("temp/benchmark",stamp,".log"), open = "wt")
  pdf(paste0("temp/benchmark",stamp,".pdf"))
  sink(zz)
  sink(zz, type = "message")
  library(doParallel)
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)

  oldwd = getwd()
  setwd(wd)

  formulas = input_formulas

  wybrane = list()
  ile_miRNA = numeric()
  for (i in 1:length(formulas)) {
    wybrane[[i]] = all.vars(formulas[[i]])
    ile_miRNA[i] = length(all.vars(formulas[[i]]))
  }

  hist(ile_miRNA, breaks = 150)
  psych::describe(ile_miRNA)

  # Wywalamy formuły z więcej niż 15 miRNA
  #formulas_old = formulas
  #formulas = formulas[which(ile_miRNA <= 15)]
  #print(formulas)

  dane = ks.wczytajmix(wd = wd, replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  if(!dir.exists("models")) { dir.create("models") }

  wyniki = data.frame(method = names(formulas))
  wyniki$SMOTE = ifelse(grepl("SMOTE", names(formulas)),"Yes","No")


  if(mxnet == T) {
    library(mxnet)
    gpu = system("nvidia-smi", intern = T)
    if (length(gpu) > 1) {
    mxctx = mx.gpu()
    print("GPU!!")
    } else {
      mxctx = mx.cpu()
      print("CPU :(")
    }


    algorytmy = c("mxnet", algorithms)
  } else {
    algorytmy = algorithms
    }

  for (ii in 1:length(algorytmy)) {
    algorytm = algorytmy[ii]
    for (i in 1:length(formulas)) {

      wyniki$miRy[i] = as.character(formulas[[i]])
      temptrain = train
      if (wyniki$SMOTE[i] == "Yes") { temptrain = train_smoted }

      # Hold-out czy CV?
      temptrainold = temptrain
      if(holdout == T) {
        fit_on = list(rs1 = 1:nrow(temptrain))
        pred_on = list(rs1 = (nrow(temptrain+1):((nrow(temptrain+1)+nrow(test)))))
        temptrain = rbind.fill(temptrain,test)
      }

      if(algorytm == "mxnet") {
          hyperparameters = expand.grid(layer1 = seq(2,14,1), layer2 = 0, layer3 = 0, activation = c('relu', 'sigmoid', 'tanh', 'softrelu'),
                                        dropout = c(0,0.05),learning.rate= c(0.00001, 0.01), momentum = c(0, 0.8, 0.99))
          train_control <- trainControl(method="cv", repeats=5, number = 10, classProbs = TRUE,verboseIter = TRUE, summaryFunction = twoClassSummary)
          if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                           classProbs = TRUE, summaryFunction = twoClassSummary) }
          model1 = caret::train(as.formula(formulas[[i]]), ctx = mxctx, optimizer = 'sgd',
                                #optimizer_params=(('learning_rate',0.1),('lr_scheduler',lr_sch)),
                                #preProc = c("center", "scale"),
                                #epoch.end.callback = mx.callback.early.stop(5,10,NULL, maximize=TRUE, verbose=TRUE),
                                #eval.data = list(data=dplyr::select(test, starts_with("hsa")),label=dplyr::select(test, Class)),
                                epoch.end.callback=mx.callback.early.stop(10, 10),
                                num.round = search_iters, data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)

      } else {
        train_control <- trainControl(method="repeatedcv", repeats=5, number = 10, search="random", classProbs = TRUE, verboseIter = TRUE,
                                      summaryFunction = twoClassSummary)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, tuneLength = search_iters)
      }

      print(model1$finalModel)
      modelname = as.numeric(Sys.time())
      print(paste0("MODELID: ", modelname))
      saveRDS(model1, paste0("models/",modelname,".RDS"))
      wyniki$modelname[i] = modelname

      t1_roc = roc(temptrainold$Class ~ predict(model1, newdata = temptrainold , type = "prob")[,2])
      wyniki[i,paste0(algorytm,"_train_ROCAUC")] = t1_roc$auc
      wyniki[i,paste0(algorytm,"_train_ROCAUC_lower95CI")] = as.character(ci(t1_roc))[1]
      wyniki[i,paste0(algorytm,"_train_ROCAUC_upper95CI")] = as.character(ci(t1_roc))[3]

      t1 = caret::confusionMatrix(predict(model1, newdata = temptrainold), as.factor(temptrainold$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_train_Accuracy")] = t1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_train_Sensitivity")] = t1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_train_Specificity")] = t1$byClass["Specificity"]

      v1 = caret::confusionMatrix(predict(model1, newdata = test), as.factor(test$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_test_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_test_Sensitivity")]  = v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_test_Specificity")]  = v1$byClass["Specificity"]

      v1 = caret::confusionMatrix(predict(model1, newdata = valid), as.factor(valid$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_valid_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_valid_Sensitivity")]= v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_valid_Specificity")] = v1$byClass["Specificity"]

      stargazer(wyniki, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)
    }
  }
  stopCluster(cl)
  write.csv(wyniki,output_file)
  setwd(oldwd)
  sink(type = "message")
  sink()
  dev.off()
  return(wyniki)
}

ks.korektaaliasowmiRNA = function(temp, species = "hsa") {
  library(foreach)
  library(doParallel)
  library(dplyr)
  miRbase_aliasy = fread("ftp://mirbase.org/pub/mirbase/CURRENT/aliases.txt.gz")
  colnames(miRbase_aliasy) = c("MIMAT","Aliasy")
  miRbase_aliasy_hsa = miRbase_aliasy %>% filter(str_detect(Aliasy, paste0(species,"*"))) %>% filter(str_detect(MIMAT, "MIMAT*"))

  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores-1) #not to overload your computer
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  temp2 = colnames(temp)
  final <- foreach(i=1:length(temp2), .combine=c) %dopar% {
    #for(i in 1:length(temp2)) {
    naz = temp2[i]
    library(data.table)
    library(stringr)
    for (ii in 1:nrow(miRbase_aliasy_hsa)) {
      temp3 = str_split(as.character(miRbase_aliasy_hsa[ii,2]), ";")
      temp4 = temp3[[1]]
      temp4 = temp4[temp4 != ""]
      if(naz %in% temp4) { naz = temp4[length(temp4)] }
    }
    naz
  }

  colnames(temp) = final
  stopCluster(cl)
  sink()
  return(temp)
}

ks.setup = function(install_pkgs = T, create_pyenv = T) {
  # apt install libcairo2-dev libopencv-dev build-essential git ninja-build ccache libopenblas-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev libtool autossh ocl-icd-opencl-dev libfreetype6-dev libssl-dev default-jdk default-jre libcurl4-openssl-dev libxml2-dev
  if(install_pkgs) {
  install.packages("BiocManager")
  BiocManager::install(c("reticulate","devtools","plyr","dplyr","edgeR","epiDisplay","rsq","MASS","Biocomb","caret","dplyr",
                       "pROC","ggplot2","DMwR", "doParallel", "Boruta", "spFSR", "varSelRF", "stringr", "psych", "C50", "randomForest",
                       "foreach","data.table", "ROSE", "deepnet", "gridExtra", "stargazer","gplots"))
  devtools::install_github("STATWORX/bounceR")
  #system("../mxnet/docs/install/install_mxnet_ubuntu_python.sh")
  }
  library(reticulate)
  library(devtools)

  # paczki nie w CRAN
  #install.packages("BiocManager")

  py_config()

  if(create_pyenv) {
  gpu = system("nvidia-smi", intern = T)
  if (length(gpu) > 1) {
    reticulate::conda_create("tensorflow",c("tensorflow-gpu","keras","mxnet-gpu"))
  } else { reticulate::conda_create("tensorflow",c("tensorflow","keras","mxnet")) }
  }

  use_condaenv("tensorflow")


}

ks.merge_formulas = function(wd = getwd(), max_miRNAs = 12) {
  oldwd = getwd()
  setwd(wd)
  formulas_files = list.files("temp","formulas*", all.files = T, full.names = T)
  formulas = list()
  formulas_names = character()
  for (i in 1:length(formulas_files)) {
    temp = readRDS(formulas_files[i])
    tempn = names(temp)
    formulas = c(formulas, temp)
    formulas_names = c(formulas_names, tempn)
  }

  temp = data.frame(name = formulas_names, formula = unlist(as.character(formulas)), stringsAsFactors = F)
  final = temp %>% dplyr::distinct()
  final$formula
  formulas = list()
  for(i in 1:nrow(final)){
    final$ile_miRNA[i] = length(all.vars(as.formula(final$formula[i])))-1
  }
  final = final %>% filter(ile_miRNA <= max_miRNAs)
  formulas_final = as.list(final$formula)
  names(formulas_final) = make.names(final$name)
  setwd(oldwd)
  saveRDS(formulas_final, "featureselection_formulas_final.RDS")
  sink()
  dev.off()
  return(formulas_final)
}

ks.benchmark.deep = function(wd = getwd(), search_iters = 200, cores = detectCores()-1, input_formulas = readRDS("featureselection_formulas_final.RDS"),
                        output_file = "benchmarkdeep.csv", algorithms = c("mxnet","mxnetAdam", "dnn"), holdout = T, stamp = as.character(as.numeric(Sys.time()))) {
  library(plyr)
  library(dplyr)
  library(edgeR)
  library(epiDisplay)
  library(rsq)
  library(MASS)
  library(Biocomb)
  library(caret)
  library(dplyr)
  library(epiDisplay)
  library(pROC)
  library(ggplot2)
  library(DMwR)
  library(stringr)
  library(psych)
  library(C50)
  library(randomForest)
  library(nnet)
  library(reticulate)
  library(deepnet)

  algorytmy = algorithms
  if(!dir.exists("temp")) { dir.create("temp") }

  use_condaenv("tensorflow")
  zz <- file(paste0("temp/benchmark",stamp,".log"), open = "wt")
  pdf(paste0("temp/benchmark",stamp,".pdf"))
  sink(zz)
  sink(zz, type = "message")
  library(doParallel)
  cl <- makePSOCKcluster(cores)
  registerDoParallel(cl)

  oldwd = getwd()
  setwd(wd)

  formulas = input_formulas

  wybrane = list()
  ile_miRNA = numeric()
  for (i in 1:length(formulas)) {
    wybrane[[i]] = all.vars(formulas[[i]])
    ile_miRNA[i] = length(all.vars(formulas[[i]]))
  }

  hist(ile_miRNA, breaks = 150)
  psych::describe(ile_miRNA)

  # Wywalamy formuły z więcej niż 15 miRNA
  #formulas_old = formulas
  #formulas = formulas[which(ile_miRNA <= 15)]
  #print(formulas)

  dane = ks.wczytajmix(wd = wd, replace_smote = F); train = dane[[1]]; test = dane[[2]]; valid = dane[[3]]; train_smoted = dane[[4]]; trainx = dane[[5]]; trainx_smoted = dane[[6]]

  if(!dir.exists("models")) { dir.create("models") }

  wyniki = data.frame(method = names(formulas))
  wyniki$SMOTE = ifelse(grepl("SMOTE", names(formulas)),"Yes","No")



  #algorytmy = c("mxnet","rf","C5.0", "nnet")
  # specyficzne train_controle
  if(sum(grepl("mxnet", algorithms)) >0) {
    library(mxnet)
    gpu = system("nvidia-smi", intern = T)
    if (length(gpu) > 1) {
      mxctx = mx.gpu() } else {
        mxctx = mx.cpu()
      }
  }

  for (ii in 1:length(algorytmy)) {
    algorytm = algorytmy[ii]
    for (i in 1:length(formulas)) {

      wyniki$miRy[i] = as.character(formulas[[i]])
      temptrain = train
      if (wyniki$SMOTE[i] == "Yes") { temptrain = train_smoted }

      temptrainold = temptrain
      if(holdout == T) {
        fit_on = list(rs1 = 1:nrow(temptrain))
        pred_on = list(rs1 = (nrow(temptrain+1):((nrow(temptrain+1)+nrow(test)))))
        temptrain = rbind.fill(temptrain,test)
      }

      if(algorytm == "mxnet") {
        hyperparameters = expand.grid(layer1 = seq(2,16,2), layer2 = seq(2,16,2), layer3 = seq(0,16,2), activation = c('relu', 'sigmoid', 'tanh', 'softrelu'),
                                      dropout = seq(0,0.1,0.05), learning.rate= c(0.00001, 0.01), momentum = c(0, 0.5, 0.9, 0.99))
        train_control <- trainControl(method="cv", repeats=5, number = 10, classProbs = TRUE,verboseIter = TRUE, summaryFunction = twoClassSummary)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary) }
        model1 = caret::train(as.formula(formulas[[i]]), ctx = mxctx, optimizer = 'sgd', num.round = search_iters, data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)
      } else if (algorytm == "mxnetAdam") {
        hyperparameters = expand.grid( layer1= seq(2,16,2), layer2 = seq(2,16,2), layer3 = seq(0,16,2), activation= c('relu', 'sigmoid', 'tanh', 'softrelu'), learningrate=c(0.00001, 0.01),
                                       beta1=0.9, beta2=0.9999, dropout=c(0,0.05,0.20) )
        train_control <- trainControl(method="cv", repeats=5, number = 10, classProbs = TRUE, verboseIter = TRUE, summaryFunction = twoClassSummary)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary) }
        model1 = caret::train(as.formula(formulas[[i]]), ctx = mxctx, optimizer = 'sgd', num.round = search_iters, data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)
      }
      else if (algorytm == "dnn") {
        hyperparameters = expand.grid( layer1 = seq(2,16,2),
                                       layer2 = seq(2,16,2), layer3 = seq(0,16,2),
                                       hidden_dropout = c(0, .1),
                                       visible_dropout = 0)
        train_control <- trainControl(method="repeatedcv", repeats=3, number = 10, classProbs = TRUE, verboseIter = TRUE, summaryFunction = twoClassSummary)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, tuneGrid = hyperparameters)
      } else
      {
        train_control <- trainControl(method="repeatedcv", repeats=3, number = 10, search="random", classProbs = TRUE,verboseIter = TRUE, summaryFunction = twoClassSummary)
        if(holdout == T) { train_control <- trainControl(method="cv", index= fit_on, indexOut = pred_on, indexFinal = fit_on[[1]], verboseIter = TRUE,
                                                         classProbs = TRUE, summaryFunction = twoClassSummary) }
        model1 = caret::train(as.formula(formulas[[i]]), data=temptrain, trControl=train_control, method=algorytm, tuneLength = search_iters)
      }

      print(model1$finalModel)
      modelname = as.numeric(Sys.time())
      print(paste0("MODELID: ", modelname))
      saveRDS(model1, paste0("models/",modelname,".RDS"))
      wyniki$modelname[i] = modelname

      t1_roc = roc(temptrainold$Class ~ predict(model1, newdata = temptrainold , type = "prob")[,2])
      wyniki[i,paste0(algorytm,"_train_ROCAUC")] = t1_roc$auc
      wyniki[i,paste0(algorytm,"_train_ROCAUC_lower95CI")] = as.character(ci(t1_roc))[1]
      wyniki[i,paste0(algorytm,"_train_ROCAUC_upper95CI")] = as.character(ci(t1_roc))[3]

      t1 = caret::confusionMatrix(predict(model1, newdata = temptrainold), as.factor(temptrainold$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_train_Accuracy")] = t1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_train_Sensitivity")] = t1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_train_Specificity")] = t1$byClass["Specificity"]

      v1 = caret::confusionMatrix(predict(model1, newdata = test), as.factor(test$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_test_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_test_Sensitivity")]  = v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_test_Specificity")]  = v1$byClass["Specificity"]

      v1 = caret::confusionMatrix(predict(model1, newdata = valid), as.factor(valid$Class), positive = "Cancer")
      wyniki[i,paste0(algorytm,"_valid_Accuracy")]  = v1$overall["Accuracy"]
      wyniki[i,paste0(algorytm,"_valid_Sensitivity")]= v1$byClass["Sensitivity"]
      wyniki[i,paste0(algorytm,"_valid_Specificity")] = v1$byClass["Specificity"]

      stargazer(wyniki, type = 'html', out = paste0("temp/wyniki",stamp,".html"), summary = F)
    }
  }
  stopCluster(cl)
  write.csv(wyniki,output_file)
  setwd(oldwd)
  sink(type = "message")
  sink()
  dev.off()
  return(wyniki)
}

ks.miRNA_signiture_overlap = function(which_formulas = c("sig","cfs"), input_formulas = readRDS("featureselection_formulas.RDS"))
{
  wybrane = list()
  for (i in 1:length(which_formulas)) {
    ktora_to = match(which_formulas[i], names(input_formulas))
    temp = as.formula(input_formulas[[ktora_to]])
    wybrane[[names(input_formulas)[ktora_to]]] = all.vars(temp)[-1]
  }
  require("gplots")
  temp = venn(wybrane)
  temp
}
