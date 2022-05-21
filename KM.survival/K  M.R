library("survival")
library("survminer")

###   需要  行  为样本，列为基因表达量，生存时间及生存事件
res.cut <- surv_cutpoint(count3, 
                         time = "time",
                         event = "event",
                         variables = "KIF20A") #可以添加多列

summary(res.cut)#查看最佳cutoff
####       cutpoint statistic
####     KIF20A    9.764      1.81
count3 <- mutate(count3,KIF20A_cutoff = ifelse(KIF20A > 9.764 ,"High" , "Low"))
fit <- survfit(Surv(time, event) ~ KIF20A_cutoff, data = count3)
ggsurvplot(fit,
           data = count3, 
           palette=c("red", "blue"), #自定义颜色
           legend.labs=c("KIF20A_High","KIF20A_Low"), #自定义标签
           risk.table = TRUE,  
           break.x.by = 6,##横坐标间隔
           pval = T) #是否展示P值
