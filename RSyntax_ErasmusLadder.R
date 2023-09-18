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
#   install.packages("tidyverse")
#   install.packages("tidyr")
#   install.packages("magrittr")
#   install.packages("gtools")
#   install.packages("ggpubr")
#   install.packages("plotrix")
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
  library(tidyverse)
  library(tidyr)
  library(magrittr)
  library(gtools)
  library(ggpubr)
  library(plotrix)
}

####ErasmusLadder Cage Enrichment project####

#set wd - working directory
rm(wd) #delete current wd
wd <- setwd("~/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/erasmusladder") #set wd
print(wd) #display wd

#read csv file, select data, convert to dataframe, recode data
#csvfile <- "erasmusladderCE.csv"
csvfile <- "erasmusladderCE.csv"

#import CSV data 
LadderRaw <- read.csv(csvfile) 
View(LadderRaw)

#convert CSV to data.frame
LadderRaw <- data.frame(LadderRaw)

#unique miceid
umid = unique(LadderRaw$MouseIdentifier);

#unique genotype/group
ugeno = unique(LadderRaw$Genotype);

#how many mice per group
setDT(LadderRaw)[, .(count = uniqueN(MouseIdentifier)), by = Genotype]

####Parameters of interest####
#Short step: step from one high rung to the next high rung
#Long step: skipping one high rung
#Jump: skipping two high rungs
#Misstep: a step forward, but the paw is placed on a low rung
#Back step: a step back from a high rung to the previous high rung 
# Value 0: means that this specific step hasn't been made during that particular trial

# Outcomes step length (correct steps):
#   2 = short high-high
#   4 = long high-high
#   6 = jump high-high

# Outcomes step length (missteps)
#   2 = short high low
#   4 = long high low
#   6 = jump high low
#----------------------------------------------------------------------------------------------------------------------------------------#

# Dataframe steps of interest: all HH (correct steps) jumps, long and short steps
{
allHH <- select(LadderRaw,Session, Trial, MouseIdentifier, Genotype, HHShortSteps, HHJumps, HHLongSteps) %>%
  mutate(totalHH = (HHShortSteps + HHJumps + HHLongSteps))
  colnames(allHH)[4] <- "Group"
  colnames(allHH) [3] <-  "AnimalID"
}

#Convert to as.factor
allHH$Session <- as.factor(allHH$Session)
allHH$Trial <- as.factor(allHH$Trial)
allHH$AnimalID <- as.factor(allHH$AnimalID)
allHH$Group <-  as.factor(allHH$Group)


### boxplot per group HH (correct) steps, normality check ###
plot1 <- ggplot(allHH,aes(x=factor(Group, level=c("Control", "CE")), y = totalHH)) +
  geom_boxplot(aes(fill = Group)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_x_discrete (name = "") +
  scale_y_continuous (name = "Correct Steps (%)") +
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


### histogram HH (correct) steps, normality check ###
{histogram1 <- ggplot(data = allHH, aes(x = totalHH)) +
    geom_histogram((aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), 
                   breaks = seq(-10, 105, by = 6)) +
    labs(title = "", x = "Correct steps (%)", y = "Frequency percentage/bin") +
    theme_bw() 
  
  plot(histogram1)}


### histogram HH (correct) steps per group, normality check ###
{histogram2 <- ggplot(data = allHH, aes(x = totalHH)) +
  geom_histogram((aes(fill = Group, y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), 
                 breaks = seq(-10, 105, by = 6)) +
    coord_cartesian(ylim = c(0, 25)) +  
  facet_wrap(factor(allHH$Group, levels=c("Control", "CE")) ~.) + 
  labs(title = "", x = "Correct steps (%)", y = "Frequency percentage/bin") +
  theme_bw() 

plot(histogram2)}


# aggregated file 1, per mouse per session
{HHaggr <- allHH %>%
    group_by(AnimalID, Group, Session) %>% 
    summarise_at(c("totalHH"), funs( 
      mean (., na.rm=T),
      median (., na.rm=T),
      #n = sum(!is.na(.)),
      #min = min (.,na.rm=T),
      #max = max (.,na.rm=T),
      sd = sd (.,na.rm=T),
      se = sd(., na.rm=T)/sqrt(sum(!is.na(.))),
      #cv = sd(., na.rm=T)/mean (., na.rm=T),
      #meanCI = ci(., na.rm=T)[1],
      #lowCI = ci(., na.rm=T)[2],
      #highCI = ci(., na.rm=T)[3],
      #sdCI = ci(., na.rm=T)[4],
      ci95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2))
  
  View(HHaggr)}


#lineplot mean correct steps all individual mice per group
{plot2 <- ggplot(HHaggr, aes(x = factor(Session), y=mean, group = AnimalID, color = Group)) +
    facet_wrap(factor(Group, levels = c("Control", "CE")) ~.) +
    geom_line() +
    scale_x_discrete (name = "Session") +
    scale_y_continuous (name = "Correct steps (%)", breaks = c(0, 40, 50, 60, 70, 80, 90)) +
    expand_limits(y=c(40,90)) +
    theme_bw() +
    ggtitle ("") +
    theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
          axis.text.y  = element_text(size=14),
          axis.line = element_line(colour = 'black', size = 0.3), 
          axis.ticks = element_line(colour = "black", size = 0.3),
          plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))
  
  plot(plot2)}


# aggregate of aggregate file 2
{HHaggr1 <- HHaggr %>%
  group_by(Group, Session) %>% 
  summarise_if(is.numeric, funs( 
              mean (., na.rm=T)))
              #median (., na.rm=T),
              #n = sum(!is.na(.)),
              #min = min (.,na.rm=T),
              #max = max (.,na.rm=T),
              #sd = sd (.,na.rm=T),
              #se = sd(., na.rm=T)/sqrt(sum(!is.na(.))),
              #cv = sd(., na.rm=T)/mean (., na.rm=T),
              #meanCI = ci(., na.rm=T)[1],
              #lowCI = ci(., na.rm=T)[2],
              #highCI = ci(., na.rm=T)[3],
              #sdCI = ci(., na.rm=T)[4],
              #ci95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2))

 View(HHaggr1)}


# lineplot, mean HH steps over 5 consecutive days
{plot3 <- ggplot(HHaggr1, aes(x=factor(Session), y=mean, group = Group, color = Group)) + 
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean - ci95,
                      ymax = mean + ci95),
                 width = 0.2,
                 size = 0.7) +
    scale_x_discrete (name = "Session") +
    scale_y_continuous (name = "Correct Steps (%)", breaks = c(0, 40, 50, 60, 70, 80, 90))  +
    expand_limits(y=c(40,90)) +
    theme_bw() +
    ggtitle ("") +
    theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
          axis.text.x  = element_text(size=14),
          axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
          axis.text.y  = element_text(size=14),
          axis.line = element_line(colour = 'black', size = 0.3), 
          axis.ticks = element_line(colour = "black", size = 0.3),
          plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

plot(plot3)}

# Linear Mixed-Effects Models correct steps ~ Group*Session
LME1<-lme(totalHH ~ Group + Session + Group*Session,
          data = allHH, 
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


#unadjusted p-values. multiple comparisons adjustment
posthoc_no <- emmeans(LME1, list(pairwise ~ Group | Session), adjust = "none") #unadjusted p-values of group per session

sumposthoc_no <- data.frame(unclass(summary(posthoc_no$`pairwise differences of Group | Session`)), check.names = FALSE, stringsAsFactors = FALSE) #created df of the posthoc_no pairwise outcome measures
sumposthoc_no


#adjusted p-values, using the Holm-Bonferroni test
posthoc_holm <- p.adjust(sumposthoc_no$p.value, method = "holm") # adjusted p-values for group per session
posthoc_holm



#p.values <- c(0.0054,0.1175,0.1206,0.0968,0.3530) # outcome measures of unadjusted p.values
#posthoc_holm <- p.adjust(p.values, method = "holm") # outcome measures adjusted p-values 
#posthoc_holm

