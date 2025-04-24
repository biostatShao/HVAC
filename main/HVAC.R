
GRS.HVAC <- function(GRS.input, traits, chr, sums_p, base_pop, target_pop,
                    base_p, base_f,target_p, target_f,cova = NULL,
                    sums_out="/data2/projects/bioinfo/zhshao/GWAS.summary/cross_pop/",
                    pheno, phe_trait, out, temp, bina = F){
  cat(sprintf("
        
        *************************************************
        *         ._  ._  .__   __   __   _____         *
        *         | |_| | | |  / /_ /  |_| .___|        *
        *         |  _  |_| | / /_ / / |_| |            *
        *         | | | |_| |/ /_ / __ |_| |___         *
        *         |_| |_| |___/  /_/ |_| |_____|        *
        * HVAC v0.0.1 (2024-12-02)                      *
        * Zhonghe Shao                                  *
        * HuaZhong University of Science and Technology *
        *************************************************
        
    \n")
  )
  
  cat(paste0("> Start HVAC method in chromesome ",chr,
             "\n> Trait is ",traits,",",
             "\n> Target population is ",target_pop,",",
             "\n> Base population is ",base_pop," \n"))
  
  out_p <- paste0(out,traits); dir.check(out_p)
  temp_p <- paste0(temp,traits); dir.check(temp_p)
  R2=T
  
  if(length(list.files(out_p,pattern = "VAC")) != 0){
    print("VAX_M1 is done")
  }else{
    if(!file.exists(paste0(sums_p,"VAX_M1"))){
      cat("--------------- Observation Layer ---------------\n")
      
      vax_methods = "VAX_M1"
      
      if(!file.exists(paste0(temp_p,"/",vax_methods,".done.",traits,"_",chr))){
        cat("Start VAX_M1 model\n")
        vax_m1 <- paste0(
          "python /data2/projects/bioinfo/zhshao/GRS/VAX/bin/vac",
          # " -l /data2/projects/bioinfo/zhshao/LD_reference/UKB/VAX/",target_pop,"/chr_",chr,
          " -l /data2/projects/bioinfo/zhshao/LD_reference/1000G/VAX/",target_pop,"/chr_",chr,
          " -s ",sums_p,traits,".fastGWA",
          " --output-dir ",temp_p,
          " --output-file-prefix ",vax_methods,".",chr,
          " --hyp-search GS",
          " --pi-steps 10",
          " --sigma-epsilon-steps 10",
          # " --lambda-min 'infer'",
          " --sumstats-format fastgwa",
          " --backend plink",
          " --validation-bfile ",target_p,target_f,chr,
          " --validation-pheno ",sums_p,traits, ".phen",
          " --validation-keep ",sums_p,traits, ".keep",
          " --temp-dir ",temp_p
          # ," > ",temp_p,"/VI",chr
        )
        system(vax_m1)
        
        prs <- paste0("/home/opt/software/plink2_linux_x86_64_20211011/plink2",
                      " --bfile ",target_p,target_f,chr,
                      " --score ",temp_p,"/",vax_methods,".",chr,".beta.gz 2 4 cols=+scoresums,-scoreavgs",
                      " --score-col-nums 6-105",
                      " --rm-dup force-first",
                      " --out ", temp_p, "/",vax_methods,"-beta-",chr,
                      " > ", temp_p, "/",vax_methods,"-beta-",chr)
        system(prs)
        
        prs <- paste0("/home/opt/software/plink2_linux_x86_64_20211011/plink2",
                      " --bfile ",target_p,target_f,chr,
                      " --score ",temp_p,"/",vax_methods,".",chr,".fit.gz 2 4 6 cols=+scoresums,-scoreavgs",
                      " --rm-dup force-first",
                      " --out ", temp_p, "/",vax_methods,"-",chr,
                      " > ", temp_p, "/",vax_methods,"-",chr)
        system(prs)
        system(paste0("touch ",temp_p,"/",vax_methods,".done.",traits,"_",chr))
      }
      
      files <- list.files(path = temp_p, pattern = paste0(vax_methods,".done.",traits,"_"))
      cat("> A total of ",length(files), " chromosome files has created\n")
      if(length(files) == 22){
        # cat("Begin merge 22 chromosome result\n")
        # options(scipen=100)
        
        rf_files = list.files(temp_p,pattern = paste0(vax_methods,"-beta")) %>% .[grep(".sscore",.)]
        res <- NULL
        for(c in 1:length(rf_files))try({
          d1 = paste0(temp_p,"/",rf_files[c]) %>% fread(.) %>%
            .[,-c(1,3,4)]
          if(c == 1){
            res <- d1
          }else{
            dd = merge(res, d1, by="IID")
            dd[,5:104] = dd[,2:101]+dd[,102:201]
            res <- dd[,1:101]
          }
        })
        fwrite(res,paste0(sums_p,vax_methods,"-beta"))
        
        rf_files = list.files(temp_p,pattern = vax_methods) %>% .[grep(".sscore",.)]
        res <- NULL
        for(c in 1:length(rf_files))try({
          d1 = paste0(temp_p,"/",rf_files[c]) %>% fread(.) %>%
            .[,c("IID", "SCORE1_SUM")]
          if(c == 1){
            res <- d1[,c("IID", "SCORE1_SUM")]
          }else{
            dd = merge(res, d1[,c("IID", "SCORE1_SUM")], by="IID")
            dd$SCORE1_SUM <- dd$SCORE1_SUM.x + dd$SCORE1_SUM.y
            res <- dd[,c("IID", "SCORE1_SUM")]
          }
        })
        fwrite(res,paste0(sums_p,vax_methods))
        
        # if(R2){
        #   phen <- fread(pheno) %>% data.frame()
        #   phen = phen[,c(colnames(phen)[1],phe_trait,cova)]
        #   data <- merge(phen, res, by.x = colnames(phen)[1], by.y = "IID") %>% data.frame()
        #   data = na.omit(data)
        #   data <- data[ceiling(nrow(data)/2):nrow(data),]
        #   
        #   if(bina == F){
        #     f1 = as.formula(paste0(phe_trait,"~",paste(cova,collapse = "+")))
        #     eta <- summary(lm(f1,data))$residuals
        #     R2 <- cor(eta,data$SCORE1_SUM)**2
        #   }else{
        #     R2 <- suppressMessages(roc(data[,phe_trait], data[,"SCORE1_SUM"]))
        #     R2 <- as.numeric(auc(R2))
        #   }
        #   fwrite(as.data.frame(1),paste0(out_p,"/",vax_methods,"-", R2))
        # }
        cat(">> Model 1 done!\n")
      }
      
    } # Model-1 VAX_target
    
  }
  
  if(length(list.files(out_p,pattern = "VAC")) != 0){
    print("VAX_M2 is done")
  }else{
    if(!file.exists(paste0(sums_p,"VAX_M2"))){
      # cat("---------------- Model 2 (VAX in base population) ----------------\n")
      
      base_traits = strsplit(traits,split = "_")[[1]][1]
      base_sums_p = gsub(pattern = traits,replacement = base_traits,x = sums_p)
      vax_methods = "VAX_M2"
      
      if(!file.exists(paste0(temp_p,"/",vax_methods,".done.",traits,"_",chr))){
        cat("Start VAX_M2 model\n")
        vax_m2 <- paste0(
          "python /data2/projects/bioinfo/zhshao/GRS/VAX/bin/vac",
          # " -l /data2/projects/bioinfo/zhshao/LD_reference/UKB/VAX/",base_pop,"/chr_",chr,
          " -l /data2/projects/bioinfo/zhshao/LD_reference/1000G/VAX/",base_pop,"/chr_",chr,
          " -s ",base_sums_p,base_traits,".fastGWA",
          " --output-dir ",temp_p,
          " --output-file-prefix ",vax_methods,".",chr,
          # " --hyp-search EM",
          " --hyp-search GS",
          " --pi-steps 10",
          " --sigma-epsilon-steps 10",
          # " --lambda-min 'infer'",
          " --sumstats-format fastgwa",
          " --backend plink",
          " --validation-bfile /data2/projects/bioinfo/zhshao/UKB_hm3eur/ukb22828_eur_hm3_chr",chr,
          " --validation-pheno ",base_sums_p,base_traits, ".phen",
          " --validation-keep ",base_sums_p,base_traits, ".keep",
          " --temp-dir ",temp_p
          # ," > ",temp_p,"/VI",chr
        )
        system(vax_m2)
        
        prs <- paste0("/home/opt/software/plink2_linux_x86_64_20211011/plink2",
                      " --bfile ",target_p,target_f,chr,
                      " --score ",temp_p,"/",vax_methods,".",chr,".beta.gz 2 4 cols=+scoresums,-scoreavgs",
                      " --score-col-nums 6-105",
                      " --rm-dup force-first",
                      " --out ", temp_p, "/",vax_methods,"-beta-",chr,
                      " > ", temp_p, "/",vax_methods,"-beta-",chr)
        system(prs)
        
        prs <- paste0("/home/opt/software/plink2_linux_x86_64_20211011/plink2",
                      " --bfile ",target_p,target_f,chr,
                      " --score ",temp_p,"/",vax_methods,".",chr,".fit.gz 2 4 6 cols=+scoresums,-scoreavgs",
                      " --rm-dup force-first",
                      " --out ", temp_p, "/",vax_methods,"-",chr,
                      " > ", temp_p, "/",vax_methods,"-",chr)
        system(prs)
        system(paste0("touch ",temp_p,"/",vax_methods,".done.",traits,"_",chr))
        
      }
      files <- list.files(path = temp_p, pattern = paste0(vax_methods,".done.",traits,"_"))
      cat("> A total of ",length(files), " chromosome files has created\n")
      if(length(files) == 22){
        # cat("Begin merge 22 chromosome result\n")
        # options(scipen=100)
        
        rf_files = list.files(temp_p,pattern = paste0(vax_methods,"-beta")) %>% .[grep(".sscore",.)]
        res <- NULL
        for(c in 1:length(rf_files))try({
          d1 = paste0(temp_p,"/",rf_files[c]) %>% fread(.) %>%
            .[,-c(1,3,4)]
          if(c == 1){
            res <- d1
          }else{
            dd = merge(res, d1, by="IID")
            dd[,5:104] = dd[,2:101]+dd[,102:201]
            res <- dd[,1:101]
          }
        })
        fwrite(res,paste0(sums_p,vax_methods,"-beta"))
        
        rf_files = list.files(temp_p,pattern = vax_methods) %>% .[grep(".sscore",.)]
        res <- NULL
        for(c in 1:length(rf_files))try({
          d1 = paste0(temp_p,"/",rf_files[c]) %>% fread(.) %>%
            .[,c("IID", "SCORE1_SUM")]
          if(c == 1){
            res <- d1[,c("IID", "SCORE1_SUM")]
          }else{
            dd = merge(res, d1[,c("IID", "SCORE1_SUM")], by="IID")
            dd$SCORE1_SUM <- dd$SCORE1_SUM.x + dd$SCORE1_SUM.y
            res <- dd[,c("IID", "SCORE1_SUM")]
          }
        })
        fwrite(res,paste0(sums_p,vax_methods))
        
        # if(R2){
        #   phen <- fread(pheno) %>% data.frame()
        #   phen = phen[,c(colnames(phen)[1],phe_trait,cova)]
        #   data <- merge(phen, res, by.x = colnames(phen)[1], by.y = "IID") %>% data.frame()
        #   data = na.omit(data)
        #   data <- data[ceiling(nrow(data)/2):nrow(data),]
        #   
        #   if(bina == F){
        #     f1 = as.formula(paste0(phe_trait,"~",paste(cova,collapse = "+")))
        #     eta <- summary(lm(f1,data))$residuals
        #     R2 <- cor(eta,data$SCORE1_SUM)**2
        #   }else{
        #     R2 <- suppressMessages(roc(data[,phe_trait], data[,"SCORE1_SUM"]))
        #     R2 <- as.numeric(auc(R2))
        #   }
        #   fwrite(as.data.frame(1),paste0(out_p,"/",vax_methods,"-", R2))
        # }
        cat(">> Model 2 done!\n")
      }
      
    } # Model-2 VAX_base
  }
  
  flen = length(list.files(path = temp_p, pattern = paste0("VAX_M1.done.",traits,"_")))
  flen = flen + length(list.files(path = temp_p, pattern = paste0("VAX_M2.done.",traits,"_")))
  
  if(length(list.files(out_p,pattern = "HVAC")) != 0){
    print("VAX_M3 is done")
  }else{
    if(flen == 44){
      if(!file.exists(paste0(sums_p,"/VAX_M3"))){
        cat("--------------- Population  Layer ---------------\n")
        
        vax_methods = "VAX_M3"
        cat("Start VAX_M3 model\n")
        res = NULL
        for(c in 1:22){
          d1 = fread(paste0(temp_p,"/VAX_M1.",c,".beta.gz"))
          d2 = fread(paste0(temp_p,"/VAX_M2.",c,".beta.gz"))
          d12 = merge(d1,d2,by = c("CHR","SNP","POS","A1","A2")) %>% as.data.frame()
          res = rbind(res,d12)
        }
        
        library(Rcpp)
        # sourceCpp("/data2/projects/bioinfo/zhshao/GRS/VAX/hvi.cpp")
        sourceCpp("/data2/projects/bioinfo/zhshao/GRS/VAX/hvi_cpp.cpp")
        
        M3_res = res[,1:5]
        for(cp in 0:99){
          mu_afr = as.numeric(res[,paste0("BETA_",cp,".x")])
          sigma_afr = as.numeric(res[,paste0("VAR_BETA_",cp,".x")])
          mu_eur = as.numeric(res[,paste0("BETA_",cp,".y")])
          sigma_eur = as.numeric(res[,paste0("VAR_BETA_",cp,".y")])
          
          # result <- vi_bayes_elbo(mu_afr, sigma_afr, mu_eur, sigma_eur, max_iter=1000, tol=1e-6)
          
          result <- vi_bayes_paper(
            mu_afr = mu_afr,
            sigma_afr = sigma_afr,
            mu_eur = mu_eur,
            sigma_eur = sigma_eur,
            max_iter = 1000,
            tol = 1e-6
          )
          M3_res = cbind(M3_res,result$beta_afr)
        }
        
        fwrite(M3_res,paste0(temp_p,"/VAX_M3.fit"),sep="\t",quote=F,row.names=F)
        
        prs <- paste0("/home/opt/software/plink2_linux_x86_64_20211011/plink2",
                      " --bfile ",target_p,"hm3_all",
                      " --score ",temp_p,"/VAX_M3.fit 2 4 cols=+scoresums,-scoreavgs",
                      " --score-col-nums 6-105",
                      " --rm-dup force-first",
                      " --out ", temp_p, "/VAX_M3",
                      " > ", temp_p, "/VAX_M3")
        system(prs)
        
        res = fread(paste0(temp_p, "/VAX_M3.sscore"))
        fwrite(res[,-c(1,3,4)],paste0(sums_p,vax_methods))
        
        # phen <- fread(pheno)
        # data <- merge(phen, res, by.x = colnames(phen)[1], by.y = "IID") %>% data.frame()
        # data = data[,which(colnames(data) %in% c(phe_trait,"SCORE1_SUM",cova))] %>% na.omit()
        # data <- data[ceiling(nrow(data)/2):nrow(data),]
        # 
        # if(bina == F){
        #   f1 = as.formula(paste0(phe_trait,"~",paste(cova,collapse = "+")))
        #   eta <- summary(lm(f1,data))$residuals
        #   R2 <- cor(eta,data$SCORE1_SUM)**2
        # }else{
        #   R2 <- suppressMessages(roc(data[,phe_trait], data[,"SCORE1_SUM"]))
        #   R2 <- as.numeric(auc(R2))
        # }
        # fwrite(as.data.frame(1),paste0(out_p,"/",vax_methods,"-", R2))
        
        cat(">> Model 3 done!\n")
      }
    } # Model-3 VAX_cross
  }
  
  # if(length(list.files(out_p,pattern = "VAC")) != 0){
  #   print("VAX_SL is done")
  # }else{
  #   
    cat("----------- Intergate models using SuperLearning ------------\n")
    M1 = fread(paste0(sums_p,"VAX_M1"))
    M2 = fread(paste0(sums_p,"VAX_M2"))
    M3 = fread(paste0(sums_p,"VAX_M3"))
    M11 = fread(paste0(sums_p,"VAX_M1-beta"))
    M21 = fread(paste0(sums_p,"VAX_M2-beta"))
    res = cbind(M1,M2[,-1],M3[,-1],M11[,-1],M21[,-1]) %>% data.table()
    colnames(res)[1] = c("FID")
    
    phen <- fread(pheno) %>% data.frame()
    phen = phen[,c(colnames(phen)[1],phe_trait,cova)]
    data <- merge(phen, res, by.x = colnames(phen)[1], by.y = "FID") %>% data.frame()
    data = na.omit(data)
    data <- data[ceiling(nrow(data)/2):nrow(data),]
    if(is.null(cova) | bina == T){
      eta <- data[,phe_trait]
    }else{
      f1 = as.formula(paste0(phe_trait,"~",paste(cova,collapse = "+")))
      eta <- summary(lm(f1,data))$residuals
    }
    
    data_train <- data[1:ceiling(nrow(data)/2),
                       -c(1:length(c(colnames(phen)[1],phe_trait,cova)))]
    eta_train <- eta[1:ceiling(nrow(data)/2)]
    data_test <- data[ceiling(nrow(data)/2):nrow(data),
                      -c(1:length(c(colnames(phen)[1],phe_trait,cova)))]
    eta_test <- eta[ceiling(nrow(data)/2):nrow(data)]
    
    if(T){
      SL.library <- c(
        "SL.bartMachine", 
        "SL.biglasso", "SL.caret", 
        "SL.caret.rpart",
        "SL.earth", "SL.gam", 
        "SL.glm", "SL.glm.interaction", 
        "SL.glmnet", "SL.ipredbagg", 
        "SL.kernelKnn", "SL.knn", "SL.ksvm", 
        "SL.lda", "SL.leekasso", "SL.lm", 
        "SL.loess", "SL.logreg", "SL.mean", 
        "SL.nnet", "SL.nnls", "SL.polymars", 
        "SL.qda", "SL.randomForest", "SL.ranger",
        "SL.ridge", "SL.rpart", "SL.rpartPrune", 
        "SL.speedglm", "SL.speedlm", "SL.step", 
        "SL.step.forward", "SL.step.interaction", 
        "SL.stepAIC", "SL.svm", "SL.template", "SL.xgboost"
      )
      
      if(bina == F){
        sl <- SuperLearner(Y = eta_train, 
                           X = data_train, 
                           family = gaussian(),
                           SL.library = c("SL.glmnet","SL.ridge"))
        pred = predict(sl,data_test,onlySL= T)
      }
      
      RR2 = NULL
      for(cc in 1:ncol(data_test)){
        if(bina == F){
          rr2 = cor(eta_test,data_test[,cc])**2
        }else{
          rr2 <- suppressMessages(roc(data[ceiling(nrow(data)/2):nrow(data),
                                           phe_trait], as.numeric(data_test[,cc])))
          rr2 = as.numeric(auc(rr2))
        }
        RR2 = c(RR2,rr2)
      }
      
      if(bina == F){
        R2 <- max(cor(eta_test,pred$pred)**2,RR2,na.rm = T)
      }else{
        # R2 <- suppressMessages(roc(data[ceiling(nrow(data)/2):nrow(data),phe_trait], as.numeric(pred$pred)))
        # R2 <- as.numeric(auc(R2))
        R2 <- max(RR2,na.rm = T)
      }
      fwrite(as.data.frame(1),paste0(out_p,"/HVAC-", R2))
    } # Super learning
    
    if(F){
      sa <- round(1/5*nrow(data_train), digits = 0)
      fold_result = NULL
      for(fold in 1:5){
        index_sa <- (1+(fold - 1)*sa):min(fold*sa, nrow(data_train))
        train = data_train[index_sa,]
        test = data_train[-index_sa,]
        Y = eta_train[index_sa]
        if(T){
          lm_coe = summary(lm(Y~.,cbind(Y,train)))$coefficients[,1]
          fold_coef = matrix(lm_coe,ncol = 1)
          
          adj_PRS = cbind(1,as.matrix(test)) %*%  fold_coef
          fold_result = rbind(fold_result, c(cor(eta_train[-index_sa],adj_PRS)**2,fold_coef))
        } # linear regression
        
        if(T){
          lambdas=exp(seq(log(0.001), log(0.1), length.out = 20))
          X = as.matrix(train)
          if(bina == F){
            lasso_model <- cv.glmnet(X,Y,alpha = 1,lambda = lambdas, nfolds = 10)
            lambda.min <- lasso_model$lambda.min
            lasso_best <- glmnet(X,Y,alpha = 1,lambda = lambda.min)
          }else{
            lasso_model <- cv.glmnet(X,Y,alpha = 1,lambda = lambdas, nfolds = 10,family = "binomial")
            lambda.min <- lasso_model$lambda.min
            lasso_best <- glmnet(X,Y,alpha = 1,lambda = lambda.min,family = "binomial")
          }
          fold_coef = coef(lasso_best) %>% matrix(.,ncol = 1)
          
          adj_PRS = cbind(1,as.matrix(test)) %*%  fold_coef
          fold_result = rbind(fold_result, c(cor(eta_train[-index_sa],adj_PRS)**2,fold_coef))
        } # LASSO
        
        if(T){
          fit <- bas.lm(Y~., data = cbind(Y,train),
                        method = "MCMC", MCMC.iterations = 10000,
                        prior = "g-prior", alpha = sqrt(nrow(cbind(Y,train))))
          sum.fit <- summary(fit)
          fold_coef = sum.fit[1:4, 1] * sum.fit[1:4, 2]
          
          adj_PRS = cbind(1,as.matrix(test)) %*%  fold_coef
          fold_result = rbind(fold_result, c(cor(eta_train[-index_sa],adj_PRS)**2,fold_coef))
        } # BMA
      }
      adj_PRS = cbind(1,as.matrix(data_test[,c("M1","M2","M3")])) %*% 
        matrix(fold_result[which.max(fold_result[,1]),-1],ncol = 1)
      
      if(bina == F){
        R2 <- cor(eta_test,adj_PRS)**2
      }else{
        R2 <- auc(roc(eta_test, adj_PRS)) %>% as.numeric()
      }
      fwrite(as.data.frame(1),paste0(out_p,"/VAX-", R2))
    }
    
    
  # }
  
  
  cat(paste0("> Done with VAC for ",traits," in chromesome ",chr,
             "\n> A total of ",flen," files have done!\n"))
}
