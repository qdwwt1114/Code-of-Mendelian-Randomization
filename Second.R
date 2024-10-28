 #设置你自己的工作路径
setwd("E:/学习/919脑功能/step1正向")

# -加载需要用到的包，没安装的需要提前安装一下
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(mr.raps)
library(ieugwasr)
library(dplyr)
library(MendelianRandomization)
library(vroom)

# 结局，把数据放在你设置的工作路径下
FileNames <-list.files(paste0(getwd()),pattern=".csv")  
exp_dat_ids <- FileNames
exps <- FileNames


# 结局
out<-vroom("E:\\学习\\919脑功能\\step1正向\\结局\\finngen_R10_ALCODELIRIUM")#读取本地结局（这里将文件名换成你自己本地结局的文件名）
out$trait <- 'ALCODELIRIUM'  #改成你所用到的结局的表型的名称
outcomeid <- out
rm(out)
head(outcomeid)
#
dir.create(path = "ALCODELIRIUM") #工作路径下新建一个文件夹存放结果

#######以下为循环代码######################
for (qaq in 1:length(exp_dat_ids)) { 
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  d3 <- try(fread(paste0(getwd(),"/",FileNames[qaq]), fill=TRUE), silent = F)
  d3 <- subset(d3, d3$pval < 5e-6)  # 暴露的p值标准
  
  if (nrow(d3) > 0) {
    d3<-as.data.frame(d3)
    d3 <- format_data(d3, type="exposure")
    
    library(ieugwasr)     
    
    d4 <- ld_clump(
      clump_kb = 500,
      clump_r2 = 0.1,    # 连锁不平衡的标准
      pop = "EUR",
      dplyr::tibble(rsid=d3$SNP, pval=d3$pval.exposure, id=d3$id.exposure),
      plink_bin = "./本地plink/plink_win64_20230116/plink.exe",
      bfile = "./本地plink/EUR/EUR"
    )

    exp_data <- subset(d3, SNP %in% d4$rsid) 
    
    if (length(exp_data[,1]) > 0) {
      outcome_dat <- merge(exp_data, outcomeid, by.x = "SNP", by.y = "SNP")
      if (length(outcome_dat[,1]) > 0) {  
        write.csv(outcome_dat, file = "d.csv")
        out_data <- read_outcome_data(
          snps = exp_data$SNP, 
          filename = "d.csv",
          sep = ","
        )
        ###排除与结局强相关
        out_data <- subset(out_data,pval.outcome>1e-5)
        
        dat <- TwoSampleMR::harmonise_data(
          exposure_dat = exp_data,
          outcome_dat = out_data)
        
        ####回文的直接去除
        dat <-subset(dat,mr_keep==TRUE)
        
        
        #####新的计算R2的代码
        dat <- get_f(dat, F_value = 10)
        
    ##############MR-PRESSO######################    
        if (nrow(dat) > 3) {
        mr_presso_res <- mr_Presso(dat, num = 1000)
        mr_presso_main <- mr_presso_pval(mr_presso_res)
        dat <- mr_presso_snp(mr_presso_res, mr_presso_main, dat, type = "data")
        #dat <-subset(dat,mr_keep==TRUE)

        resMRPRESSO=mr_presso_res[["Main MR results"]]
        resMRPRESSO
        global_test_p <- mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
        se1=sqrt(((resMRPRESSO[1,3])^2)/qchisq(Pval_raw <- resMRPRESSO[1,6],1,lower.tail=F))
        se2=sqrt(((beta_cor <- resMRPRESSO[2,3])^2)/qchisq(Pval_cor <- resMRPRESSO[2,6],1,lower.tail=F))
        resMRPRESSO <- resMRPRESSO %>%
          dplyr::mutate(se = c(se1,se2))
        
        #生成globe test_P和outliers
        outliers <- dat$SNP[mr_presso_res[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]]
        outliers = as.data.frame(outliers)
        outliers <- paste(outliers$outliers, collapse = ", ")
        global_test_p = as.data.frame(global_test_p)
        resMRPRESSO
        TTT <- as.data.frame(resMRPRESSO)
        TTT
        openxlsx::write.xlsx(TTT, paste0("ALCODELIRIUM/",exp,"-MR-PRESSO.xlsx"))
        } 
        dat <-subset(dat,mr_keep==TRUE)
        res <-choose_MR(dat=dat)
        res_hete <-purrr::map(.x=seq_along(res),
                              .f=~res[[.x]][[1]])
        res_hete <-do.call(rbind,res_hete)
        res_hete
        res1 <- generate_odds_ratios(res[[1]][[2]])
        res1
        res1$estimate <- paste0(
          format(round(res1$or, 2), nsmall = 2), " (", 
          format(round(res1$or_lci95, 2), nsmall = 2), "-",
          format(round(res1$or_uci95, 2), nsmall = 2), ")")
        res1
        print(paste0(exp,"_SNP数_",res1$nsnp[1]))
        resdata <- dat
        openxlsx::write.xlsx(dat,file = paste0("ALCODELIRIUM/",exp,"-dat.xlsx"), rowNames = FALSE)
        
        openxlsx::write.xlsx(res1,paste0("ALCODELIRIUM/",exp,"-res.xlsx"))
        
        
        ######steiger检验######
        res_steiger <- steiger_test(dat) 
        
        
        # ##二分类结局power计算
        N <- dat$samplesize.outcome[1]
        alpha <- 0.05
        R2xz <- sum(dat$R2)
        K <- (dat$ncase.outcome[1] / dat$ncontrol.outcome[1])
        if (length(dat[, 1]) == 1) {
          OR <- res1 %>% filter(method == "Wald ratio") %>% pull(or)
        } else if (length(dat[, 1]) > 1) {
          OR <-
            res1 %>%filter(grepl("Inverse variance weighted", method)) %>%pull(or)
        }
        epower = NA
        power <- results_binary(N, alpha, R2xz, K, OR, epower)
        
        library(magrittr)
        # Main result 
        res3 <- res1#[1:5,]
        res3 <- res3[,-c(10:14)]
        res4 <- tidyr::pivot_wider(
          res3,names_from ="method",names_vary = "slowest",
          values_from = c("b","se","pval","estimate") )
        
        col_names <- colnames(res4)
        
        new_col_names <- gsub("\\(.*\\)", "", col_names)
        
        colnames(res4) <- new_col_names
        
        ##steiger
        res_steiger2 <- dplyr::select(res_steiger,
                                      correct_causal_direction,steiger_pval)
        #Power
        power2 <- tidyr::pivot_wider(
          power,names_from ="Parameter",names_vary = "slowest",
          values_from = c("Value","Description") )
        power2 <- power2[,1]
        
        
        # Merge
        res_ALL <- cbind(res4, res_steiger2,power2
        )
        res_ALL$F <- mean(dat$F,na.rm=TRUE)
        res_ALL$R2 <- sum(dat$R2)
        
        
        
        if (length(dat[, 1]) > 0 && length(dat[, 1]) <= 2) {
          write.csv(res_ALL, file = paste0("ALCODELIRIUM/", exp, "1.csv"), row.names = FALSE)
        } else {
          # res_hete <- TwoSampleMR::mr_heterogeneity(dat)
          res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
          res_leaveone <- mr_leaveoneout(dat)
          
          p1 <- mr_scatter_plot(res[[1]][[2]], dat)
          pdf(paste0("ALCODELIRIUM/", exp, "_scatter.pdf"))
          print(p1[[1]])
          dev.off()
          
          res_single <- mr_singlesnp(dat,all_method)
          p2 <- mr_forest_plot(res_single)
          pdf(paste0("ALCODELIRIUM/", exp, "_forest.pdf"))
          print(p2[[1]])
          dev.off()
          
          p3 <- mr_funnel_plot(res_single)
          pdf(paste0("ALCODELIRIUM/", exp, "_funnel.pdf"))
          print(p3[[1]])
          dev.off()
          
          res_loo <- mr_leaveoneout(dat)
          pdf(paste0("ALCODELIRIUM/", exp, "_leave_one_out.pdf"))
          print(mr_leaveoneout_plot(res_loo))
          dev.off()
          
          
          
          res_hete2 <- tidyr::pivot_wider(
            res_hete, names_from = "method", names_vary = "slowest",
            values_from = c("Q", "Q_df", "Q_pval")
          ) %>% 
            dplyr::select(-id.exposure, -id.outcome, -outcome, -exposure)
          res_hete2 <- res_hete2[, 4:6]
          
          res_plei2 <- dplyr::select(res_plei, egger_intercept, se, pval)
          
          res_ALL <- cbind(res_ALL, res_hete2, res_plei2)
          
          if (nrow(dat) > 3) {
  
          res_ALL$outliers <- outliers
          res_ALL <- cbind(res_ALL, global_test_p)
          }
          write.csv(res_ALL, file = paste0("ALCODELIRIUM/", exp, ".csv"), row.names = FALSE)
        }
      }
    }
  }
}



#最后一步才运行，把导出的路径改成自己的工作路径  

library(dplyr)

# 获取所有CSV文件的路径
csv_files <- list.files("./ALCODELIRIUM", 
                        pattern = "\\.csv$", full.names = TRUE)      

# 读取第一个文件
combined_df <- read.csv(csv_files[1], stringsAsFactors = FALSE)

# 检查并修正数据结构
if (!is.null(combined_df$exposure)) {      
  combined_df$exposure <- as.character(combined_df$exposure)      
}   
if (!is.null(combined_df$global_test_p)) {      
  combined_df$global_test_p <- as.character(combined_df$global_test_p)      
}  

# 遍历并合并剩余文件
for (i in 2:length(csv_files)) {      
  temp_df <- read.csv(csv_files[i], stringsAsFactors = FALSE)      
  
  # 转换列为字符型，确保一致性
  if (!is.null(temp_df$exposure)) {      
    temp_df$exposure <- as.character(temp_df$exposure)      
  }    
  if (!is.null(temp_df$global_test_p)) {      
    temp_df$global_test_p <- as.character(temp_df$global_test_p)      
    temp_df$global_test_p[temp_df$global_test_p == "<"] <- ""      
  }      
  
  # 确保所有列的数据类型一致
  temp_df <- temp_df %>% mutate(across(everything(), as.character))
  combined_df <- combined_df %>% mutate(across(everything(), as.character))
  
  # 合并数据框
  combined_df <- bind_rows(combined_df, temp_df)      
}

  #计算Fdr后的P（Benjamini-Hochberg方法）
combined_df <- combined_df %>%
  mutate(fdr_Inverse.variance.weighted. = p.adjust(pval_Inverse.variance.weighted., method = "fdr")) %>%
  relocate(fdr_Inverse.variance.weighted., .after = pval_Inverse.variance.weighted.)
write.csv(combined_df,"ALCODELIRIUMres.csv")






