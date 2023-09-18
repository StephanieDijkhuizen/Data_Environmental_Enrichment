####install libraries####
# { 
#   install.packages("nlme")
#   install.packages("data.table")
#   install.packages("emmeans")
#   install.packages("lmerTest")
#   install.packages("lattice")
#   install.packages("plyr")
#   install.packages("ggplot2")
#   install.packages("Renext")
#   install.packages("ggsn")
#   install.packages("doBy")
#   install.packages("Rmisc")
#   install.packages("dplyr")
#   install.packages("gmodels")
#   install.packages("googlesheets")
#   install.packages("httpuv")
#   install.packages("RCurl")
#   install.packages("gsheet")
#   install.packages("standardize")
#   install.packages("lme4")
#   install.packages("Matrix")
#   }

####used libraries####
{
  library(nlme)
  library(data.table)
  library(emmeans)
  library(lmerTest)
  library(lattice)
  library(plyr)
  library(ggplot2)
  library(Renext)
  library(ggsn)
  library(doBy)
  library(Rmisc)
  library(dplyr)
  library(gmodels)
  library(googlesheets)
  library(httpuv)
  library(RCurl)
  library(gsheet)
  library(standardize)
  library(lme4)
  library(Matrix)
}

#Dropbox (Neurasmus BV) Grip Strength Test

# set wd - working directory
rm(wd)
wd <- setwd("~/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/Rotarod")
wd <- setwd("~/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/Rotarod")
print(wd)


# read csv file, select data, convert to dataframe, recode data
#csvfile <- "RR Total Data.csv"
csvfile <- "RR_Total_Data.csv"
RotaRaw <- read.csv(csvfile)
View(RotaRaw)

write.csv(RotaRaw, file = "RotaRaw.csv")

#RotaRaw excl Day5, due to the difference in RPM
RotaRaw <- RotaRaw %>% 
  filter(!(Day == 5))
View(RotaRaw)

write.csv(RotaRaw, file = "RotaRaw_no_day5.csv")

# convert CSV to data.frame
RotaRaw <- data.frame(RotaRaw)

#Convert to factor
RotaRaw$Day <- as.factor(RotaRaw$Day)
RotaRaw$Block <- as.factor(RotaRaw$Block)
RotaRaw$AnimalID <- as.factor(RotaRaw$AnimalID)


# plot1, geom_boxplot + geom_dotplot ALL data points of Day 1 - Day 4
plot1 <- ggplot(RotaRaw, aes(x=factor(Genotype, level=c("Control", "CE")),y=Latencytofall)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1.1, dotsize=0.6, width = 3, aes(fill = Genotype)) +
  scale_x_discrete (name = "Genotype") +
  scale_y_continuous (name = "Latency to fall (s.)") +
  theme_bw() +
  ggtitle ("") +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.y  = element_text(size=14),
        axis.line = element_line(colour = 'black', size = 0.3), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

plot(plot1)

# Linear Mixed-Effects Models
LME1<-lme(Latencytofall ~ Genotype + Day + Genotype*Day,
          data = RotaRaw, 
          correlation = NULL,
          random = ~1 | AnimalID, 
          method = "REML", 
          na.action = na.exclude)
summary(LME1)
anova(LME1)
VarCorr(LME1)
plot(LME1)
qqnorm(LME1)

sum <- (summary(LME1))
an <- (anova(LME1))
vc <- (VarCorr(LME1))


# unadjusted p-values. multiple comparisons adjustment
posthoc_no <- emmeans(LME1, list(pairwise ~ Genotype | Day), adjust = "none") #unadjusted p-values of group per session

sumposthoc_no <- data.frame(unclass(summary(posthoc_no$`pairwise differences of Genotype | Day` )), check.names = FALSE, stringsAsFactors = FALSE) #created df of the posthoc_no pairwise outcome measures
sumposthoc_no

#adjusted p-values using Holm-bonferroni 
posthoc_holm <- p.adjust(sumposthoc_no$p.value, method = "holm") # adjusted p-values for group per session
posthoc_holm


# unadjusted p-values. multiple comparisons adjustment
#posthoc_no <- emmeans(LME1, list(pairwise ~ Genotype | Day), adjust = "none)")
#posthoc_no
#p.no <- (posthoc_no)

# adjusted p-values. multiple comparisons adjustment
#p.values <- c(0.0003,0.0026,0.0411,0.0934) # outcome measures of unadjusted p.values
#p.adjust(p.values, method = "holm")


# make aggregate file 1 which is aggregate
rm(RotaRawAggr)
{RotaRawAggr <- RotaRaw %>% 
    group_by (Genotype, AnimalID, Day, Sex) %>% 
    summarise_at(c("Latencytofall"), funs( 
  mean (., na.rm=T),
  median (., na.rm=T),
  #n = sum(!is.na(.)),
  #min = min (.,na.rm=T),
  #max = max (.,na.rm=T),
  sd = sd (.,na.rm=T),
  se = sd(., na.rm=T)/sqrt(sum(!is.na(.))), 
  #cv = sd(., na.rm=T)/mean (., na.rm=T),
  # meanCI = ci (., na.rm=T)[1],
  #lowCI = ci (., na.rm=T)[2],
  #highCI = ci (., na.rm=T)[3],
  # sdCI = ci (., na.rm=T)[4]
  ci95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2))
}

View(RotaRawAggr)

write.csv(RotaRawAggr, file = "RotaRawAggr.csv")

# plot2, geom_boxplot + geom_dotplot
plot2 <- ggplot(RotaRawAggr, aes(x=factor(Genotype, level=c("Control", "CE")),y=mean)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1.1, dotsize=0.8, width = 3, aes(fill = Genotype)) +
  scale_x_discrete (name = "") +
  scale_y_continuous (name = "Latency to fall (s.)") +
  theme_bw() +
  ggtitle ("") +
  theme(legend.text = element_text(size = 12)) +
  guides(col = guide_legend(ncol = 1, byrow = TRUE, keywidth = 0.1, keyheight = 0.1, default.unit = "inch")) +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.y  = element_text(size=14),
        axis.line = element_line(colour = 'black', size = 0.3), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

plot(plot2)


# make aggregate1 file 2
rm(RotaRawAggr1)
{RotaRawAggr1 <- RotaRaw %>% 
    group_by (Genotype, Day, Block) %>% 
    summarise_at (c("Latencytofall"), funs( 
  mean (., na.rm=T),
  median (., na.rm=T),
  #n = sum(!is.na(.)),
  #min = min (.,na.rm=T),
  #max = max (.,na.rm=T),
  sd = sd (.,na.rm=T),
  se = sd(., na.rm=T)/sqrt(sum(!is.na(.))), 
  #cv = sd(., na.rm=T)/mean (., na.rm=T),
  # meanCI = ci (., na.rm=T)[1],
  #lowCI = ci (., na.rm=T)[2],
  #highCI = ci (., na.rm=T)[3],
  # sdCI = ci (., na.rm=T)[4]
  ci95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2))
}

View(RotaRawAggr1)

write.csv(RotaRawAggr1, file = "RotaRawAggr1.csv")


#plot3, line plot mean of each block per day
 plot3 <- ggplot(RotaRawAggr1, aes(x=factor(Block),y=mean,group = Genotype,color = Genotype)) + 
  facet_grid(~ Day) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = mean - ci95,
                    ymax = mean + ci95),
                    width = 0.2,
                    size = 0.7) +
  scale_x_discrete (name = "Blocks/Day") +
  scale_y_continuous (name = "Latency to fall (s.)") +
  expand_limits(y=c(0,300)) +
  theme_bw() +
  ggtitle ("Rotarod Test") +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
       axis.text.x  = element_text(size=14),
       axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
       axis.text.y  = element_text(size=14),
       axis.line = element_line(colour = 'black', size = 0.3), 
       axis.ticks = element_line(colour = "black", size = 0.3),
       plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

  plot(plot3)

 
  

# make aggregate1 file 3 which is aggregate of aggregate
{RotaRawAggr2 <- RotaRawAggr1 %>% 
    group_by (Genotype, Day) %>% 
    summarise_if (is.numeric, funs( 
      mean (., na.rm=T)))
      #median (., na.rm=T),
      #n = sum(!is.na(.)),
      #min = min (.,na.rm=T),
      #max = max (.,na.rm=T),
      #sd = sd (.,na.rm=T),
      #se = sd(., na.rm=T)/sqrt(sum(!is.na(.))), 
      #cv = sd(., na.rm=T)/mean (., na.rm=T),
      # meanCI = ci (., na.rm=T)[1],
      #lowCI = ci (., na.rm=T)[2],
      #highCI = ci (., na.rm=T)[3],
      # sdCI = ci (., na.rm=T)[4]
      #i95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2))}
  
  View(RotaRawAggr2)}
  
  write.csv(RotaRawAggr2, file = "RotaRawAggr2.csv")
  

#plot4, line plot mean per day
  plot4 <- ggplot(RotaRawAggr2, aes(x=factor(Day),y=mean,group = Genotype,color = Genotype)) + 
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - ci95,
                      ymax = mean + ci95),
                  width = 0.2,
                  size = 0.7) +
    scale_x_discrete (name = "Day") +
    scale_y_continuous (name = "Latency to fall (s.)") +
    expand_limits(y=c(0,300)) +
    theme_bw() +
    ggtitle ("Rotarod Test") +
    theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
          axis.text.y  = element_text(size=14),
          axis.line = element_line(colour = 'black', size = 0.3), 
          axis.ticks = element_line(colour = "black", size = 0.3),
          plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))
  
  plot(plot4)
  


  #LME2 per session in CE group
  colnames(RotaRaw)
  unique(RotaRaw$Genotype)
  
  
  #subset CE group
  RotaRaw_CE <- subset (RotaRaw, Genotype == "CE")
  unique(RotaRaw_CE$Genotype)
  
  
  #LME CE group
  LME2<-lme(Latencytofall ~ Day,
            data = RotaRaw_CE, 
            correlation = NULL,
            random = ~1 | AnimalID, 
            method = "REML", 
            na.action = na.exclude)
  summary(LME2)
  anova(LME2)
  VarCorr(LME2)
  plot(LME2)
  qqnorm(LME2)
  
  
  #unadjusted p-values. multiple comparisons adjustment
  posthoc_no <- emmeans(LME2, list(pairwise ~ Day), adjust = "none") #unadjusted p-values of group per session
  
  sumposthoc_no <- data.frame(unclass(summary(posthoc_no$`pairwise differences of Day`)), check.names = FALSE, stringsAsFactors = FALSE) #created df of the posthoc_no pairwise outcome measures
  sumposthoc_no
  
  
  #adjusted p-values using Holm-bonferroni 
  posthoc_holm <- p.adjust(sumposthoc_no$p.value, method = "holm") # adjusted p-values for group per session
  posthoc_holm
  
  #LME3 per session in Control group
  colnames(RotaRaw)
  unique(RotaRaw$Genotype)
  
  
  #subset control group
  RotaRaw_Control <- subset (RotaRaw, Genotype == "Control")
  unique(RotaRaw_Control$Genotype)
  
  
  #LME control group
  LME3<-lme(Latencytofall ~ Day,
            data = RotaRaw_Control, 
            correlation = NULL,
            random = ~1 | AnimalID, 
            method = "REML", 
            na.action = na.exclude)
  summary(LME3)
  anova(LME3)
  VarCorr(LME3)
  plot(LME3)
  qqnorm(LME3)
  
  
  #unadjusted p-values. multiple comparisons adjustment
  posthoc_no <- emmeans(LME3, list(pairwise ~ Day), adjust = "none") #unadjusted p-values of group per session
  
  sumposthoc_no <- data.frame(unclass(summary(posthoc_no$`pairwise differences of Day`)), check.names = FALSE, stringsAsFactors = FALSE) #created df of the posthoc_no pairwise outcome measures
  sumposthoc_no
  
  #adjusted p-values using Holm-bonferroni 
  posthoc_holm <- p.adjust(sumposthoc_no$p.value, method = "holm") # adjusted p-values for group per session
  posthoc_holm
  
  