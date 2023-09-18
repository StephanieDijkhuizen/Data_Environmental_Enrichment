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
  library (gmodels)
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
wd <- setwd("~/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/Griptest")
wd <- setwd("~/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/Griptest")
print(wd)


# read csv file, select data, convert to dataframe, recode data
#csvfile <- "Griptest Total Data.csv"
csvfile <- "Griptest Total Data.csv"

GripRaw <- read.csv(csvfile)
View(GripRaw)


# convert CSV to data.frame
GripRaw <- data.frame(GripRaw)

write.csv(GripRaw, file = "GripRaw.csv")

#plot1
plot1 <- ggplot(GripRaw, aes(x=factor(Genotype, level=c("Control", "CE")),y=Strength)) +
  geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=1.1, dotsize=0.7, width = 3, aes(fill = Genotype)) +
  scale_x_discrete (name = "Genotype") +
  scale_y_continuous (name = "Strength (N)") +
  expand_limits(y=c(0,200)) +
  theme_bw() +
  ggtitle ("Grip Strength Test") +
  theme(legend.text = element_text(size = 12)) +
  guides(col = guide_legend(ncol = 1, byrow = TRUE, keywidth = 0.1, keyheight = 0.1, default.unit = "inch")) +
  theme(axis.title.x = element_text(size=16, margin = margin(t = 5, r = 0, b = 0, l = 0)),
        axis.text.x  = element_text(size=14),
        axis.title.y = element_text(size=16, margin = margin(t = 0, r = 10, b = 0, l = 0)), 
        axis.text.y  = element_text(size=14),
        axis.line = element_line(colour = 'black', size = 0.3), 
        axis.ticks = element_line(colour = "black", size = 0.3),
        plot.title = element_text(size = 16, face = "bold", margin = margin(t = 0, r = 0, b = 15, l = 0)))

plot(plot1)

LME1<-lme(Strength~Genotype,
          data = GripRaw, 
          correlation = NULL,
          random = ~1 | Animal.ID, 
          method = "REML", 
          na.action = na.exclude)

summary(LME1)
anova(LME1)
VarCorr(LME1)
plot(LME1)
qqnorm(LME1)



