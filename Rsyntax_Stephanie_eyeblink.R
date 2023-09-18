#######################################################################
##      PREPARE DATASET AND CREATE FOLDER IN WORKING DIRECTORY      ###
#######################################################################

  
# used libraries
{
  library(nlme)
  library(data.table)
  library(emmeans)
  library(lmerTest)
  library(emmeans)
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
}

#### data analysis cage enrichment acquisition ####

# set wd - working directory
rm(wd)
wd <- setwd("/Users/stephaniedijkhuizen1/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/eyeblink/ISI 250 ms/") # fill out the working directory, where your data is stored
wd <- setwd("/Users/stephaniedijkhuizen1/Dropbox (Neurasmus BV)/Projects/Cage Enrichment/analysis/eyeblink/ISI 500 ms/timeX750ms") # fill out the working directory, where your data is stored
print(wd)

# read csv file, select data, convert to dataframe, recode data
csvfile <- "outcomes_trial_by_trial_c57bl6cage_ce_AND_c57bl6cage_controle.csv" # fill out the the name of the trial by trial outcome file, which was exported from the mblink program
ebcRaw <- read.csv(csvfile)
nrow(ebcRaw)

# subset of ebcRaw df, only considered valid trials
rm(ebcSelect)
ebcSelect <- subset (ebcRaw, validity == "valid") # & (is.na(CRonset) # | CRonset < 750) & (is.na(CRpeaktime) | CRpeaktime > 600))
nrow(ebcSelect)

# subset CR.onset 
{ebcSelect <- subset (ebcSelect, session_nr >7)
nrow(ebcSelect)
colnames(ebcSelect)
}

# histogram CR.onset
{histogramCR.onset <- ggplot(ebcSelect, aes(x = CR.onset, color=genotype, fill=genotype)) +
  geom_histogram((aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), alpha=0.5, binwidth=(13) ) +  
  facet_grid(genotype~.) +
  scale_x_continuous(name = "CR onset", limits = (c(0,1500))) +  
  scale_y_continuous(name="Frequency percentage/bin", limits = (c(0,20))) + 
  theme_classic ()
histogramCR.onset
}

# subset CR.peaktime
{ebcSelect <- subset (ebcSelect, session_nr >7)
nrow(ebcSelect)
colnames(ebcSelect)

histogramCR.peaktime <- ggplot(ebcSelect, aes(x = CR.peaktime, color=genotype, fill=genotype)) +
  geom_histogram((aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), alpha=0.5, binwidth=(13)) + 
  facet_grid(genotype~.) +
  scale_x_continuous(name= "CR peaktime", limits = (c(0,1500))) +
  scale_y_continuous(name="Frequency percentage/bin", limits = (c(0,15))) + 
  theme_classic ()
print(histogramCR.peaktime)


lme1 <- lme (CR.peaktime ~ genotype,
                  data = ebcSelect, 
                  random = ~ 1 | mouse_id,
                  method = "REML", 
                  na.action=na.exclude) 

summary (lme1)
anova (lme1)
qqnorm(lme1)
}

# subset perfect timed CR's
{ebcSelect1 <- subset (ebcSelect, CR.peaktime > 750) # select only trials that do have a CR | >500ms for 250ms ISI and >750 for 500ms ISI
ebcSelect1$CRperc_perfect <- 0 # set all to zero
ebcSelect1$CRperc_perfect[ebcSelect1$CR.peaktime > 900 & ebcSelect1$CR.peaktime < 1100] <- 100 # set perfect to 100
nrow(ebcSelect1)
colnames(ebcSelect1)
}

# make aggregate file 1 which is aggregate
Aggr <- ebcSelect1 %>% 
    group_by (genotype) %>% 
    summarise_at(c("CRperc_perfect"), funs( 
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

View(Aggr)

# boxplot CR.peaktime CR's >700 <800 for 250ms ISI | >900 <1100 for 500ms ISI
boxplot_CR.peaktime <- ggplot(data=ebcSelect1, aes(x=genotype, y=CR.peaktime)) +
  geom_boxplot(aes(color=genotype)) +
  theme_classic () +
  stat_boxplot(geom ='errorbar', width = 0.3, aes(color=genotype)) 
print(boxplot_CR.peaktime)

# boxplot perfect timed CR's
boxplot_CRperc_perfect <- ggplot(data=aggr, aes(x=genotype, y=CRperc_perfect)) +
  geom_boxplot(aes(color=genotype)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5) +
  theme_classic() +
  scale_y_continuous(name= "perfect CR", limits = (c(0,100))) +
  stat_boxplot(geom = 'errorbar', width = 0.3, aes(color=genotype))
print(boxplot_CRperc_perfect)

#violine plot perfect timed CR's
violinplot_CRperfect <- ggplot(aggr, aes(x=genotype, y=CRperc_perfect)) +
  geom_violin(aes(fill = genotype)) 
violinplot_CRperfect

# LME2 CR.peaktime
lme2 <- lme (CR.peaktime ~ genotype,
             data = ebcSelect1,
             random = ~ 1 | mouse_id,
             method = "REML",
             na.action=na.exclude)

summary (lme2)
anova (lme2)
qqnorm(lme2)

# LME3 CRperc_perfect
lme3 <- lme (CRperc_perfect ~ genotype,
        data = ebcSelect1,
        random = ~ 1 | mouse_id,
        method = "REML",
        na.action=na.exclude)

summary (lme3)
anova (lme3)
qqnorm(lme3)


# convert CSV to data.frame
rm(ebcAll)
ebcAll <- data.frame(ebcSelect)
View(ebcAll)
nrow(ebcAll)

# check column names in df
colnames (ebcAll)

# unique
unique(ebcAll$genotype)
unique(ebcAll$mouse_id)
unique(ebcAll$trialtype)

setDT(ebcAll)[, .(count = uniqueN(mouse_id)), by = genotype] # count animals per genotype

# rename variables
{
ebcAll$URnorm <- ebcAll$UR.norm
ebcAll$UR.norm <- NULL

ebcAll$CRonset <- ebcAll$CR.onset
ebcAll$CR.onset <- NULL

ebcAll$CRpeakamp <- ebcAll$CR.peakamp
ebcAll$CR.peakamp <- NULL

ebcAll$CRpeaktime <- ebcAll$CR.peaktime
ebcAll$CR.peaktime <- NULL

ebcAll$CRhpamp <- ebcAll$CR.hp.amp
ebcAll$CR.hp.amp <- NULL

ebcAll$CRhptime <- ebcAll$CR.hp.time
ebcAll$CR.hp.time <- NULL

ebcAll$ISIpeakamp <- ebcAll$ISI.peakamp
ebcAll$ISI.peakamp <- NULL

ebcAll$ISIpeaktime <- ebcAll$ISI.peaktime
ebcAll$ISI.peaktime <- NULL

ebcAll$ISIhpamp <- ebcAll$ISI.hp.amp
ebcAll$ISI.hp.amp <- NULL

ebcAll$ISIhptime <- ebcAll$ISI.hp.time
ebcAll$ISI.hp.time <- NULL

ebcAll$AmpattX<- ebcAll$Amp.at.time.X
ebcAll$Amp.at.time.X <- NULL
}



# calculate a few new variables
{
ebcAll$log2AmpattX <- log2(ebcAll$AmpattX^2)
ebcAll$log10AmpattX <- log10(ebcAll$AmpattX^2)
ebcAll$sqrtAmpattX <- sqrt(ebcAll$AmpattX^2)

ebcAll$log2ISIhpamp <- log2(ebcAll$ISIhpamp^2)
ebcAll$log10ISIhpamp <- log10(ebcAll$ISIhpamp^2)
ebcAll$sqrtISIhpamp <- sqrt(ebcAll$ISIhpamp^2)

ebcAll$log2ISIpeakamp <- log2(ebcAll$ISIpeakamp^2)
ebcAll$log10ISIpeakamp <- log10(ebcAll$ISIpeakamp^2)
ebcAll$sqrtISIpeakamp <- sqrt(ebcAll$ISIpeakamp^2)

ebcAll$log2CRhpamp <- log2(ebcAll$CRhpamp^2)
ebcAll$log10CRhpamp <- log10(ebcAll$CRhpamp^2)
ebcAll$sqrtCRhpamp <- sqrt(ebcAll$CRhpamp^2)

ebcAll$log2CRpeakamp <- log2(ebcAll$CRpeakamp^2)
ebcAll$log10CRpeakamp <- log10(ebcAll$CRpeakamp^2)
ebcAll$sqrtCRpeakamp <- sqrt(ebcAll$CRpeakamp^2)

ebcAll$CRpeaktime_abs  <- ((abs(ebcAll$CRpeaktime-(750))/250)*-10)+10
ebcAll$CRpeakamp_abs   <- (ebcAll$CRpeakamp/150)*10
ebcAll$CRpeakscore     <- (ebcAll$CRpeaktime_abs)*(ebcAll$CRpeakamp_abs)
  
ebcAll$CRhptime_abs  <- ((abs(ebcAll$CRhptime-(750))/250)*-10)+10
ebcAll$CRhpamp_abs   <- (ebcAll$CRhpamp/150)*10
ebcAll$CRhpscore     <- (ebcAll$CRhptime_abs)*(ebcAll$CRhpamp_abs)
  
ebcAll$ISIpeaktime_abs  <- ((abs(ebcAll$ISIpeaktime-(750))/250)*-10)+10
ebcAll$ISIpeakamp_abs   <- (ebcAll$ISIpeakamp/150)*10
ebcAll$ISIpeakscore     <- (ebcAll$ISIpeaktime_abs)*(ebcAll$ISIpeakamp_abs)
  
ebcAll$ISIhptime_abs  <- ((abs(ebcAll$ISIhptime-(750))/250)*-10)+10
ebcAll$ISIhpamp_abs   <- (ebcAll$ISIhpamp/150)*10
ebcAll$ISIhpscore     <- (ebcAll$ISIhptime_abs)*(ebcAll$ISIhpamp_abs)

ebcAll$CRperc5.peak[ebcAll$ISIpeakamp < .05] <- 0
ebcAll$CRperc5.peak[ebcAll$ISIpeakamp >= .05] <- 100

ebcAll$CRperc10.peak[ebcAll$ISIpeakamp < .10] <- 0
ebcAll$CRperc10.peak[ebcAll$ISIpeakamp >= .10] <- 100

ebcAll$CRperc15.peak[ebcAll$ISIpeakamp < .15] <- 0
ebcAll$CRperc15.peak[ebcAll$ISIpeakamp >= .15] <- 100

ebcAll$CRperc20.peak[ebcAll$ISIpeakamp < .20] <- 0
ebcAll$CRperc20.peak[ebcAll$ISIpeakamp >= .20] <- 100

ebcAll$CRperc5.hp[ebcAll$ISIhpamp < .05] <- 0
ebcAll$CRperc5.hp[ebcAll$ISIhpamp >= .05] <- 100

ebcAll$CRperc10.hp[ebcAll$ISIhpamp < .10] <- 0
ebcAll$CRperc10.hp[ebcAll$ISIhpamp >= .10] <- 100

ebcAll$CRperc15.hp[ebcAll$ISIhpamp < .15] <- 0
ebcAll$CRperc15.hp[ebcAll$ISIhpamp >= .15] <- 100

ebcAll$CRperc20.hp[ebcAll$ISIhpamp < .20] <- 0
ebcAll$CRperc20.hp[ebcAll$ISIhpamp >= .20] <- 100

ebcAll$CRperc5.timeX[ebcAll$ISIhpamp < .05] <- 0
ebcAll$CRperc5.timeX[ebcAll$ISIhpamp >= .05] <- 100

ebcAll$CRperc10.timeX[ebcAll$ISIhpamp < .10] <- 0
ebcAll$CRperc10.timeX[ebcAll$ISIhpamp >= .10] <- 100

ebcAll$CRperc15.timeX[ebcAll$ISIhpamp < .15] <- 0
ebcAll$CRperc15.timeX[ebcAll$ISIhpamp >= .15] <- 100

ebcAll$CRperc20.timeX[ebcAll$ISIhpamp < .20] <- 0
ebcAll$CRperc20.timeX[ebcAll$ISIhpamp >= .20] <- 100

ebcAll$ISIperc_t600_t900.hp[ebcAll$ISIhptime < 600 | ebcAll$ISIhptime > 900] <- 0
ebcAll$ISIperc_t600_t900.hp[ebcAll$ISIhptime >= 600 & ebcAll$ISIhptime <= 900] <- 100

ebcAll$ISIperc_t650_t850.hp[ebcAll$ISIhptime < 650 | ebcAll$ISIhptime > 850] <- 0
ebcAll$ISIperc_t650_t850.hp[ebcAll$ISIhptime >= 650 & ebcAll$ISIhptime <= 850] <- 100

ebcAll$ISIperc_t700_t800.hp[ebcAll$ISIhptime < 700 | ebcAll$ISIhptime > 800] <- 0
ebcAll$ISIperc_t700_t800.hp[ebcAll$ISIhptime >= 700 & ebcAll$ISIhptime <= 800] <- 100

ebcAll$ISIperc_t725_t775.hp[ebcAll$ISIhptime < 725 | ebcAll$ISIhptime > 775] <- 0
ebcAll$ISIperc_t725_t775.hp[ebcAll$ISIhptime >= 725 & ebcAll$ISIhptime <= 775] <- 100

ebcAll$ISIperc_t650_t1000.hp[ebcAll$ISIhptime < 650 | ebcAll$ISIhptime > 1000] <- 0
ebcAll$ISIperc_t650_t1000.hp[ebcAll$ISIhptime >= 650 & ebcAll$ISIhptime <= 1000] <- 100

ebcAll$ISIperc_t600_t900.peak[ebcAll$ISIpeaktime < 600 | ebcAll$ISIpeaktime > 900] <- 0
ebcAll$ISIperc_t600_t900.peak[ebcAll$ISIpeaktime >= 600 & ebcAll$ISIpeaktime <= 900] <- 100

ebcAll$ISIperc_t650_t850.peak[ebcAll$ISIpeaktime < 650 | ebcAll$ISIpeaktime > 850] <- 0
ebcAll$ISIperc_t650_t850.peak[ebcAll$ISIpeaktime >= 650 & ebcAll$ISIpeaktime <= 850] <- 100

ebcAll$ISIperc_t700_t800.peak[ebcAll$ISIpeaktime < 700 | ebcAll$ISIpeaktime > 800] <- 0
ebcAll$ISIperc_t700_t800.peak[ebcAll$ISIpeaktime >= 700 & ebcAll$ISIpeaktime <= 800] <- 100

ebcAll$ISIperc_t725_t775.peak[ebcAll$ISIpeaktime < 725 | ebcAll$ISIpeaktime > 775] <- 0
ebcAll$ISIperc_t725_t775.peak[ebcAll$ISIpeaktime >= 725 & ebcAll$ISIpeaktime <= 775] <- 100

ebcAll$ISIperc_t650_t1000.peak[ebcAll$ISIpeaktime < 650 | ebcAll$ISIpeaktime > 1000] <- 0
ebcAll$ISIperc_t650_t1000.peak[ebcAll$ISIpeaktime >= 650 & ebcAll$ISIpeaktime <= 1000] <- 100

ebcAll$CRperc_t600_t900.hp[ebcAll$ISIhptime < 600 | ebcAll$ISIhptime > 900 | ebcAll$ISIhpamp < .10] <- 0
ebcAll$CRperc_t600_t900.hp[ebcAll$ISIhptime >= 600 & ebcAll$ISIhptime <= 900 & ebcAll$ISIhpamp >= .10] <- 100

ebcAll$CRperc_t650_t850.hp[ebcAll$ISIhptime < 650 | ebcAll$ISIhptime > 850 | ebcAll$ISIhpamp < .10] <- 0
ebcAll$CRperc_t650_t850.hp[ebcAll$ISIhptime >= 650 & ebcAll$ISIhptime <= 850 & ebcAll$ISIhpamp >= .10] <- 100

ebcAll$CRperc_t700_t800.hp[ebcAll$ISIhptime < 700 | ebcAll$ISIhptime > 800 | ebcAll$ISIhpamp < .10] <- 0
ebcAll$CRperc_t700_t800.hp[ebcAll$ISIhptime >= 700 & ebcAll$ISIhptime <= 800 & ebcAll$ISIhpamp >= .10] <- 100

ebcAll$CRperc_t725_t775.hp[ebcAll$ISIhptime < 725 | ebcAll$ISIhptime > 775 | ebcAll$ISIhpamp < .10] <- 0
ebcAll$CRperc_t725_t775.hp[ebcAll$ISIhptime >= 725 & ebcAll$ISIhptime <= 775 & ebcAll$ISIhpamp >= .10] <- 100

ebcAll$CRperc_t650_t1000.hp[ebcAll$ISIhptime < 650 | ebcAll$ISIhptime > 1000 | ebcAll$ISIhpamp < .10] <- 0
ebcAll$CRperc_t650_t1000.hp[ebcAll$ISIhptime >= 650 & ebcAll$ISIhptime <= 1000 & ebcAll$ISIhpamp >= .10] <- 100

ebcAll$CRperc_t600_t900.peak[ebcAll$ISIpeaktime < 600 | ebcAll$ISIpeaktime > 900 | ebcAll$ISIpeakamp < .10] <- 0
ebcAll$CRperc_t600_t900.peak[ebcAll$ISIpeaktime >= 600 & ebcAll$ISIpeaktime <= 900 & ebcAll$ISIpeakamp >= .10] <- 100

ebcAll$CRperc_t650_t850.peak[ebcAll$ISIpeaktime < 650 | ebcAll$ISIpeaktime > 850 | ebcAll$ISIpeakamp < .10] <- 0
ebcAll$CRperc_t650_t850.peak[ebcAll$ISIpeaktime >= 650 & ebcAll$ISIpeaktime <= 850 & ebcAll$ISIpeakamp >= .10] <- 100

ebcAll$CRperc_t700_t800.peak[ebcAll$ISIpeaktime < 700 | ebcAll$ISIpeaktime > 800 | ebcAll$ISIpeakamp < .10] <- 0
ebcAll$CRperc_t700_t800.peak[ebcAll$ISIpeaktime >= 700 & ebcAll$ISIpeaktime <= 800 & ebcAll$ISIpeakamp >= .10] <- 100

ebcAll$CRperc_t725_t775.peak[ebcAll$ISIpeaktime < 725 | ebcAll$ISIpeaktime > 775 | ebcAll$ISIpeakamp < .10] <- 0
ebcAll$CRperc_t725_t775.peak[ebcAll$ISIpeaktime >= 725 & ebcAll$ISIpeaktime <= 775 & ebcAll$ISIpeakamp >= .10] <- 100

ebcAll$CRperc_t650_t1000.peak[ebcAll$ISIpeaktime < 650 | ebcAll$ISIpeaktime > 1000 | ebcAll$ISIpeakamp < .10] <- 0
ebcAll$CRperc_t650_t1000.peak[ebcAll$ISIpeaktime >= 650 & ebcAll$ISIpeaktime <= 1000 & ebcAll$ISIpeakamp >= .10] <- 100

# recode into factors, i.e. categorical values
ebcAll$construct <- as.factor(ebcAll$construct)
ebcAll$genotype <- as.factor(ebcAll$genotype)
ebcAll$mouse_id <- as.factor(ebcAll$mouse_id)
ebcAll$session_nr <- as.factor(ebcAll$session_nr)
ebcAll$block_nr <- as.factor(ebcAll$block_nr)
ebcAll$trial_nr <- as.factor(ebcAll$trial_nr)
ebcAll$trialtype <- as.factor(ebcAll$trialtype)
ebcAll$validity <- as.factor(ebcAll$validity)

# define amplitude covariate
ebcAll$amplitude_covariate <- as.numeric(ebcAll$AmpattX)

# recode session into numeric for models 1, 2, and 3.
ebcAll$session_nr2 <- as.numeric(ebcAll$session_nr)
ebcAll$block_nr2 <- as.numeric(ebcAll$block_nr)
}


# make aggretate file 
rm(ebcAggr)
ebcAggr <- {ddply (ebcAll, .(genotype, mouse_id, session_nr), summarize, 
                   CRperc5.hp_N = sum (!is.na(CRperc5.hp)),CRperc5.hp_mean = mean (CRperc5.hp, na.rm=TRUE), CRperc5.hp_median = median (CRperc5.hp, na.rm=TRUE),
                   CRperc5.hp_sd = sd(CRperc5.hp, na.rm=TRUE), CRperc5.hp_se = CRperc5.hp_sd/sqrt(CRperc5.hp_N), CRperc5.hp_cv = CRperc5.hp_sd/CRperc5.hp_mean,
                   CRperc5.hp_CV2 = CV2 ((as.vector(CRperc5.hp))[!is.na(as.vector(CRperc5.hp))]),
                   
                   CRperc10.hp_N = sum (!is.na(CRperc10.hp)),CRperc10.hp_mean = mean (CRperc10.hp, na.rm=TRUE), CRperc10.hp_median = median (CRperc10.hp, na.rm=TRUE),
                   CRperc10.hp_sd = sd(CRperc10.hp, na.rm=TRUE), CRperc10.hp_se = CRperc10.hp_sd/sqrt(CRperc10.hp_N), CRperc10.hp_cv = CRperc10.hp_sd/CRperc10.hp_mean, 
                   CRperc10.hp_CV2 = CV2 ((as.vector(CRperc10.hp))[!is.na(as.vector(CRperc10.hp))]),
                   
                   CRperc15.hp_N = sum (!is.na(CRperc15.hp)),CRperc15.hp_mean = mean (CRperc15.hp, na.rm=TRUE), CRperc15.hp_median = median (CRperc15.hp, na.rm=TRUE),
                   CRperc15.hp_sd = sd(CRperc15.hp, na.rm=TRUE), CRperc15.hp_se = CRperc15.hp_sd/sqrt(CRperc15.hp_N), CRperc15.hp_cv = CRperc15.hp_sd/CRperc15.hp_mean, 
                   CRperc15.hp_CV2 = CV2 ((as.vector(CRperc15.hp))[!is.na(as.vector(CRperc15.hp))]),
                   
                   CRperc20.hp_N = sum (!is.na(CRperc20.hp)),CRperc20.hp_mean = mean (CRperc20.hp, na.rm=TRUE), CRperc20.hp_median = median (CRperc20.hp, na.rm=TRUE),
                   CRperc20.hp_sd = sd(CRperc20.hp, na.rm=TRUE), CRperc20.hp_se = CRperc20.hp_sd/sqrt(CRperc20.hp_N), CRperc20.hp_cv = CRperc20.hp_sd/CRperc20.hp_mean, 
                   CRperc20.hp_CV2 = CV2 ((as.vector(CRperc20.hp))[!is.na(as.vector(CRperc20.hp))]),
                   
                   CRperc5.peak_N = sum (!is.na(CRperc5.peak)),CRperc5.peak_mean = mean (CRperc5.peak, na.rm=TRUE), CRperc5.peak_median = median (CRperc5.peak, na.rm=TRUE),
                   CRperc5.peak_sd = sd(CRperc5.peak, na.rm=TRUE), CRperc5.peak_se = CRperc5.peak_sd/sqrt(CRperc5.peak_N), CRperc5.peak_cv = CRperc5.peak_sd/CRperc5.peak_mean, 
                   CRperc5.peak_CV2 = CV2 ((as.vector(CRperc5.peak))[!is.na(as.vector(CRperc5.peak))]),
                   
                   CRperc10.peak_N = sum (!is.na(CRperc10.peak)),CRperc10.peak_mean = mean (CRperc10.peak, na.rm=TRUE), CRperc10.peak_median = median (CRperc10.peak, na.rm=TRUE),
                   CRperc10.peak_sd = sd(CRperc10.peak, na.rm=TRUE), CRperc10.peak_se = CRperc10.peak_sd/sqrt(CRperc10.peak_N), CRperc10.peak_cv = CRperc10.peak_sd/CRperc10.peak_mean, 
                   CRperc10.peak_CV2 = CV2 ((as.vector(CRperc10.peak))[!is.na(as.vector( CRperc10.peak))]),
                   
                   CRperc15.peak_N = sum (!is.na(CRperc15.peak)),CRperc15.peak_mean = mean (CRperc15.peak, na.rm=TRUE), CRperc15.peak_median = median (CRperc15.peak, na.rm=TRUE),
                   CRperc15.peak_sd = sd(CRperc15.peak, na.rm=TRUE), CRperc15.peak_se = CRperc15.peak_sd/sqrt(CRperc15.peak_N), CRperc15.peak_cv = CRperc15.peak_sd/CRperc15.peak_mean, 
                   CRperc15.peak_CV2 = CV2 ((as.vector(CRperc15.peak))[!is.na(as.vector(CRperc15.peak))]),
                   
                   CRperc20.peak_N = sum (!is.na(CRperc20.peak)),CRperc20.peak_mean = mean (CRperc20.peak, na.rm=TRUE), CRperc20.peak_median = median (CRperc20.peak, na.rm=TRUE),
                   CRperc20.peak_sd = sd(CRperc20.peak, na.rm=TRUE), CRperc20.peak_se = CRperc20.peak_sd/sqrt(CRperc20.peak_N), CRperc20.peak_cv = CRperc20.peak_sd/CRperc20.peak_mean,
                   CRperc20.peak_CV2 = CV2 ((as.vector(CRperc20.peak  ))[!is.na(as.vector(CRperc20.peak  ))]),
                   
                   CRonset_N = sum (!is.na(CRonset)),CRonset_mean = mean (CRonset, na.rm=TRUE), CRonset_median = median (CRonset, na.rm=TRUE),
                   CRonset_sd = sd(CRonset, na.rm=TRUE), CRonset_se = CRonset_sd/sqrt(CRonset_N), CRonset_cv = CRonset_sd/CRonset_mean, 
                   CRonset_CV2 = CV2 ((as.vector( CRonset ))[!is.na(as.vector( CRonset ))]),
                   
                   CRhpamp_N = sum (!is.na(CRhpamp)),CRhpamp_mean = mean (CRhpamp, na.rm=TRUE), CRhpamp_median = median (CRhpamp, na.rm=TRUE),
                   CRhpamp_sd = sd(CRhpamp, na.rm=TRUE), CRhpamp_se = CRhpamp_sd/sqrt(CRhpamp_N), 
                   CRhpamp_cv = sqrt((CRhpamp_sd*CRhpamp_sd)/(CRhpamp_mean*CRhpamp_mean)), 
                   CRhpamp_CV2 =CV2 ((as.vector( CRhpamp ))[!is.na(as.vector( CRhpamp ))]),
                   
                   CRhptime_N = sum (!is.na(CRhptime)),CRhptime_mean = mean (CRhptime, na.rm=TRUE), CRhptime_median = median (CRhptime, na.rm=TRUE),
                   CRhptime_sd = sd(CRhptime, na.rm=TRUE), CRhptime_se = CRhptime_sd/sqrt(CRhptime_N), CRhptime_cv = CRhptime_sd/CRhptime_mean, 
                   CRhptime_CV2 = CV2 ((as.vector( CRhptime ))[!is.na(as.vector( CRhptime ))]),
                   
                   CRpeakamp_N = sum (!is.na(CRpeakamp)),CRpeakamp_mean = mean (CRpeakamp, na.rm=TRUE), CRpeakamp_median = median (CRpeakamp, na.rm=TRUE),
                   CRpeakamp_sd = sd(CRpeakamp, na.rm=TRUE), CRpeakamp_se = CRpeakamp_sd/sqrt(CRpeakamp_N), 
                   CRpeakamp_cv = sqrt((CRpeakamp_sd*CRpeakamp_sd)/(CRpeakamp_mean*CRpeakamp_mean)), 
                   CRpeakamp_CV2 = CV2 ((as.vector(CRpeakamp))[!is.na(as.vector(CRpeakamp))]),
                   
                   CRpeaktime_N = sum (!is.na(CRpeaktime)),CRpeaktime_mean = mean (CRpeaktime, na.rm=TRUE), CRpeaktime_median = median (CRpeaktime, na.rm=TRUE),
                   CRpeaktime_sd = sd(CRpeaktime, na.rm=TRUE), CRpeaktime_se = CRpeaktime_sd/sqrt(CRpeaktime_N), CRpeaktime_cv = CRpeaktime_sd/CRpeaktime_mean, 
                   CRpeaktime.hp_CV2 = CV2 ((as.vector( CRpeaktime ))[!is.na(as.vector(CRpeaktime  ))]),
                   
                   ISIhpamp_N = sum (!is.na(ISIhpamp)),ISIhpamp_mean = mean (ISIhpamp, na.rm=TRUE), ISIhpamp_median = median (ISIhpamp, na.rm=TRUE),
                   ISIhpamp_sd = sd(ISIhpamp, na.rm=TRUE), ISIhpamp_se = ISIhpamp_sd/sqrt(ISIhpamp_N), 
                   ISIhpamp_cv = sqrt((ISIhpamp_sd*ISIhpamp_sd)/(ISIhpamp_mean*ISIhpamp_mean)), 
                   ISIhpamp_CV2 = CV2 ((as.vector( ISIhpamp  ))[!is.na(as.vector(  ISIhpamp ))]),
                   
                   ISIhptime_N = sum (!is.na(ISIhptime)),ISIhptime_mean = mean (ISIhptime, na.rm=TRUE), ISIhptime_median = median (ISIhptime, na.rm=TRUE),
                   ISIhptime_sd = sd(ISIhptime, na.rm=TRUE), ISIhptime_se = ISIhptime_sd/sqrt(ISIhptime_N), ISIhptime_cv = ISIhptime_sd/ISIhptime_mean, 
                   ISIhptime_CV2 = CV2 ((as.vector( ISIhptime ))[!is.na(as.vector(ISIhptime  ))]),
                   
                   ISIpeakamp_N = sum (!is.na(ISIpeakamp)),ISIpeakamp_mean = mean (ISIpeakamp, na.rm=TRUE), ISIpeakamp_median = median (ISIpeakamp, na.rm=TRUE),
                   ISIpeakamp_sd = sd(ISIpeakamp, na.rm=TRUE), ISIpeakamp_se = ISIpeakamp_sd/sqrt(ISIpeakamp_N), 
                   ISIpeakamp_cv = sqrt((ISIpeakamp_sd*ISIpeakamp_sd)/(ISIpeakamp_mean*ISIpeakamp_mean)), 
                   ISIpeakamp_CV2 = CV2 ((as.vector( ISIpeakamp ))[!is.na(as.vector( ISIpeakamp ))]),
                   
                   ISIpeaktime_N = sum (!is.na(ISIpeaktime)),ISIpeaktime_mean = mean (ISIpeaktime, na.rm=TRUE), ISIpeaktime_median = median (ISIpeaktime, na.rm=TRUE),
                   ISIpeaktime_sd = sd(ISIpeaktime, na.rm=TRUE), ISIpeaktime_se = ISIpeaktime_sd/sqrt(ISIpeaktime_N), ISIpeaktime_cv = ISIpeaktime_sd/ISIpeaktime_mean, 
                   ISIpeaktime_CV2 = CV2 ((as.vector(ISIpeaktime ))[!is.na(as.vector( ISIpeaktime ))]),
                   
                   AmpattX_N = sum (!is.na(AmpattX)),AmpattX_mean = mean (AmpattX, na.rm=TRUE), AmpattX_median = median (AmpattX, na.rm=TRUE),
                   AmpattX_sd = sd(AmpattX, na.rm=TRUE), AmpattX_se = AmpattX_sd/sqrt(AmpattX_N), 
                   AmpattX_cv = sqrt((AmpattX_sd*AmpattX_sd)/(AmpattX_mean*AmpattX_mean)),
                   AmpattX_CV2 = CV2 ((as.vector(AmpattX  ))[!is.na(as.vector(AmpattX  ))]),
                   
                   CRpeakamp_abs_N = sum (!is.na(CRpeakamp_abs)),CRpeakamp_abs_mean = mean (CRpeakamp_abs, na.rm=TRUE), CRpeakamp_abs_median = median (CRpeakamp_abs, na.rm=TRUE),
                   CRpeakamp_abs_sd = sd(CRpeakamp_abs, na.rm=TRUE), CRpeakamp_abs_se = CRpeakamp_abs_sd/sqrt(CRpeakamp_abs_N), 
                   CRpeakamp_abs_cv = sqrt((CRpeakamp_abs_sd*CRpeakamp_abs_sd)/(CRpeakamp_abs_mean*CRpeakamp_abs_mean)), 
                   CRpeakamp_abs_CV2 = CV2 ((as.vector(CRpeaktime_abs  ))[!is.na(as.vector( CRpeaktime_abs ))]),
                   
                   CRpeaktime_abs_N = sum (!is.na(CRpeaktime_abs)),CRpeaktime_abs_mean = mean (CRpeaktime_abs, na.rm=TRUE), CRpeaktime_abs_median = median (CRpeaktime_abs, na.rm=TRUE),
                   CRpeaktime_abs_sd = sd(CRpeaktime_abs, na.rm=TRUE), CRpeaktime_abs_se = CRpeaktime_abs_sd/sqrt(CRpeaktime_abs_N), CRpeaktime_abs_cv = CRpeaktime_abs_sd/CRpeaktime_abs_mean, 
                   CRpeaktime_abs_CV2 = CV2 ((as.vector(CRpeaktime_abs  ))[!is.na(as.vector(CRpeaktime_abs  ))]),
                   
                   CRpeakscore_N = sum (!is.na(CRpeakscore)),CRpeakscore_mean = mean (CRpeakscore, na.rm=TRUE), CRpeakscore_median = median (CRpeakscore, na.rm=TRUE),
                   CRpeakscore_sd = sd(CRpeakscore, na.rm=TRUE), CRpeakscore_se = CRpeakscore_sd/sqrt(CRpeakscore_N), CRpeakscore_cv = CRpeakscore_sd/CRpeakscore_mean, 
                   CRpeakscore_CV2 = CV2 ((as.vector(CRpeakscore  ))[!is.na(as.vector(CRpeakscore  ))]),
                   
                   CRhpamp_abs_N = sum (!is.na(CRhpamp_abs)),CRhpamp_abs_mean = mean (CRhpamp_abs, na.rm=TRUE), CRhpamp_abs_median = median (CRhpamp_abs, na.rm=TRUE),
                   CRhpamp_abs_sd = sd(CRhpamp_abs, na.rm=TRUE), CRhpamp_abs_se = CRhpamp_abs_sd/sqrt(CRhpamp_abs_N), 
                   CRhpamp_abs_cv = sqrt((CRhpamp_abs_sd*CRhpamp_abs_sd)/(CRhpamp_abs_mean*CRhpamp_abs_mean)), 
                   CRhpamp_abs_CV2 = CV2 ((as.vector( CRhpamp_abs ))[!is.na(as.vector( CRhpamp_abs ))]),
                   
                   CRhptime_abs_N = sum (!is.na(CRhptime_abs)),CRhptime_abs_mean = mean (CRhptime_abs, na.rm=TRUE), CRhptime_abs_median = median (CRhptime_abs, na.rm=TRUE),
                   CRhptime_abs_sd = sd(CRhptime_abs, na.rm=TRUE), CRhptime_abs_se = CRhptime_abs_sd/sqrt(CRhptime_abs_N), CRhptime_abs_cv = CRhptime_abs_sd/CRhptime_abs_mean, 
                   CRhptime_abs_CV2 = CV2 ((as.vector(CRhptime_abs  ))[!is.na(as.vector( CRhptime_abs ))]),
                   
                   CRhpscore_N = sum (!is.na(CRhpscore)),CRhpscore_mean = mean (CRhpscore, na.rm=TRUE), CRhpscore_median = median (CRhpscore, na.rm=TRUE),
                   CRhpscore_sd = sd(CRhpscore, na.rm=TRUE), CRhpscore_se = CRhpscore_sd/sqrt(CRhpscore_N), CRhpscore_cv = CRhpscore_sd/CRhpscore_mean,
                   CRhpscore_CV2 = CV2 ((as.vector( CRhpscore ))[!is.na(as.vector( CRhpscore ))]),
                   
                   ISIpeakamp_abs_N = sum (!is.na(ISIpeakamp_abs)),ISIpeakamp_abs_mean = mean (ISIpeakamp_abs, na.rm=TRUE), ISIpeakamp_abs_median = median (ISIpeakamp_abs, na.rm=TRUE),
                   ISIpeakamp_abs_sd = sd(ISIpeakamp_abs, na.rm=TRUE), ISIpeakamp_abs_se = ISIpeakamp_abs_sd/sqrt(ISIpeakamp_abs_N), 
                   ISIpeakamp_abs_cv = sqrt((ISIpeakamp_abs_sd*ISIpeakamp_abs_sd)/(ISIpeakamp_abs_mean*ISIpeakamp_abs_mean)), 
                   ISIpeakamp_abs_CV2 = CV2 ((as.vector( ISIpeakamp_abs ))[!is.na(as.vector( ISIpeakamp_abs ))]),
                   
                   ISIpeaktime_abs_N = sum (!is.na(ISIpeaktime_abs)),ISIpeaktime_abs_mean = mean (ISIpeaktime_abs, na.rm=TRUE), ISIpeaktime_abs_median = median (ISIpeaktime_abs, na.rm=TRUE),
                   ISIpeaktime_abs_sd = sd(ISIpeaktime_abs, na.rm=TRUE), ISIpeaktime_abs_se = ISIpeaktime_abs_sd/sqrt(ISIpeaktime_abs_N), ISIpeaktime_abs_cv = ISIpeaktime_abs_sd/ISIpeaktime_abs_mean, 
                   ISIpeaktime_abs_CV2 = CV2 ((as.vector(ISIpeaktime_abs  ))[!is.na(as.vector(ISIpeaktime_abs  ))]),
                   
                   ISIpeakscore_N = sum (!is.na(ISIpeakscore)),ISIpeakscore_mean = mean (ISIpeakscore, na.rm=TRUE), ISIpeakscore_median = median (ISIpeakscore, na.rm=TRUE),
                   ISIpeakscore_sd = sd(ISIpeakscore, na.rm=TRUE), ISIpeakscore_se = ISIpeakscore_sd/sqrt(ISIpeakscore_N), ISIpeakscore_cv = ISIpeakscore_sd/ISIpeakscore_mean, 
                   ISIpeakscore_CV2 = CV2 ((as.vector(  ISIpeakscore ))[!is.na(as.vector(  ISIpeakscore ))]),
                   
                   ISIhpamp_abs_N = sum (!is.na(ISIhpamp_abs)),ISIhpamp_abs_mean = mean (ISIhpamp_abs, na.rm=TRUE), ISIhpamp_abs_median = median (ISIhpamp_abs, na.rm=TRUE),
                   ISIhpamp_abs_sd = sd(ISIhpamp_abs, na.rm=TRUE), ISIhpamp_abs_se = ISIhpamp_abs_sd/sqrt(ISIhpamp_abs_N), 
                   ISIhpamp_abs_cv = sqrt((ISIhpamp_abs_sd*ISIhpamp_abs_sd)/(ISIhpamp_abs_mean*ISIhpamp_abs_mean)), 
                   ISIhpamp_abs_CV2 = CV2 ((as.vector( ISIhpamp_abs ))[!is.na(as.vector( ISIhpamp_abs ))]),
                   
                   ISIhptime_abs_N = sum (!is.na(ISIhptime_abs)),ISIhptime_abs_mean = mean (ISIhptime_abs, na.rm=TRUE), ISIhptime_abs_median = median (ISIhptime_abs, na.rm=TRUE),
                   ISIhptime_abs_sd = sd(ISIhptime_abs, na.rm=TRUE), ISIhptime_abs_se = ISIhptime_abs_sd/sqrt(ISIhptime_abs_N), ISIhptime_abs_cv = ISIhptime_abs_sd/ISIhptime_abs_mean, 
                   ISIhptime_abs_CV2 = CV2 ((as.vector(  ISIhptime_abs ))[!is.na(as.vector( ISIhptime_abs  ))]),
                   
                   ISIhpscore_N = sum (!is.na(ISIhpscore)),ISIhpscore_mean = mean (ISIhpscore, na.rm=TRUE), ISIhpscore_median = median (ISIhpscore, na.rm=TRUE),
                   ISIhpscore_sd = sd(ISIhpscore, na.rm=TRUE), ISIhpscore_se = ISIhpscore_sd/sqrt(ISIhpscore_N), ISIhpscore_cv = ISIhpscore_sd/ISIhpscore_mean, 
                   ISIhpscore_CV2 = CV2 ((as.vector(ISIhpscore  ))[!is.na(as.vector(ISIhpscore  ))]),
                   
                   ISIperc_t600_t900.hp_N = sum (!is.na(ISIperc_t600_t900.hp)),ISIperc_t600_t900.hp_mean = mean (ISIperc_t600_t900.hp, na.rm=TRUE), ISIperc_t600_t900.hp_median = median (ISIperc_t600_t900.hp, na.rm=TRUE),
                   ISIperc_t600_t900.hp_sd = sd(ISIperc_t600_t900.hp, na.rm=TRUE), ISIperc_t600_t900.hp_se = ISIperc_t600_t900.hp_sd/sqrt(ISIperc_t600_t900.hp_N), ISIperc_t600_t900.hp_cv = ISIperc_t600_t900.hp_sd/ISIperc_t600_t900.hp_mean, 
                   ISIperc_t600_t900.hp_CV2 = CV2 ((as.vector(  ISIperc_t600_t900.hp ))[!is.na(as.vector(  ISIperc_t600_t900.hp ))]),
                   
                   ISIperc_t650_t850.hp_N = sum (!is.na(ISIperc_t650_t850.hp)),ISIperc_t650_t850.hp_mean = mean (ISIperc_t650_t850.hp, na.rm=TRUE), ISIperc_t650_t850.hp_median = median (ISIperc_t650_t850.hp, na.rm=TRUE),
                   ISIperc_t650_t850.hp_sd = sd(ISIperc_t650_t850.hp, na.rm=TRUE), ISIperc_t650_t850.hp_se = ISIperc_t650_t850.hp_sd/sqrt(ISIperc_t650_t850.hp_N), ISIperc_t650_t850.hp_cv = ISIperc_t650_t850.hp_sd/ISIperc_t650_t850.hp_mean, 
                   ISIperc_t650_t850.hp_CV2 = CV2 ((as.vector(ISIperc_t650_t850.hp  ))[!is.na(as.vector(ISIperc_t650_t850.hp  ))]),
                   
                   ISIperc_t700_t800.hp_N = sum (!is.na(ISIperc_t700_t800.hp)),ISIperc_t700_t800.hp_mean = mean (ISIperc_t700_t800.hp, na.rm=TRUE), ISIperc_t700_t800.hp_median = median (ISIperc_t700_t800.hp, na.rm=TRUE),
                   ISIperc_t700_t800.hp_sd = sd(ISIperc_t700_t800.hp, na.rm=TRUE), ISIperc_t700_t800.hp_se = ISIperc_t700_t800.hp_sd/sqrt(ISIperc_t700_t800.hp_N), ISIperc_t700_t800.hp_cv = ISIperc_t700_t800.hp_sd/ISIperc_t700_t800.hp_mean, 
                   ISIperc_t700_t800.hp_CV2 = CV2 ((as.vector(ISIperc_t700_t800.hp  ))[!is.na(as.vector( ISIperc_t700_t800.hp ))]),
                   
                   ISIperc_t725_t775.hp_N = sum (!is.na(ISIperc_t725_t775.hp)),ISIperc_t725_t775.hp_mean = mean (ISIperc_t725_t775.hp, na.rm=TRUE), ISIperc_t725_t775.hp_median = median (ISIperc_t725_t775.hp, na.rm=TRUE),
                   ISIperc_t725_t775.hp_sd = sd(ISIperc_t725_t775.hp, na.rm=TRUE), ISIperc_t725_t775.hp_se = ISIperc_t725_t775.hp_sd/sqrt(ISIperc_t725_t775.hp_N), ISIperc_t725_t775.hp_cv = ISIperc_t725_t775.hp_sd/ISIperc_t725_t775.hp_mean, 
                   ISIperc_t725_t775.hp_CV2 = CV2 ((as.vector(  ISIperc_t725_t775.hp ))[!is.na(as.vector(  ISIperc_t725_t775.hp ))]),
                   
                   ISIperc_t650_t1000.hp_N = sum (!is.na(ISIperc_t650_t1000.hp)),ISIperc_t650_t1000.hp_mean = mean (ISIperc_t650_t1000.hp, na.rm=TRUE), ISIperc_t650_t1000.hp_median = median (ISIperc_t650_t1000.hp, na.rm=TRUE),
                   ISIperc_t650_t1000.hp_sd = sd(ISIperc_t650_t1000.hp, na.rm=TRUE), ISIperc_t650_t1000.hp_se = ISIperc_t650_t1000.hp_sd/sqrt(ISIperc_t650_t1000.hp_N), ISIperc_t650_t1000.hp_cv = ISIperc_t650_t1000.hp_sd/ISIperc_t650_t1000.hp_mean, 
                   ISIperc_t650_t1000.hp_CV2 = CV2 ((as.vector( ISIperc_t650_t1000.hp ))[!is.na(as.vector( ISIperc_t650_t1000.hp ))]),
                   
                   ISIperc_t600_t900.peak_N = sum (!is.na(ISIperc_t600_t900.peak)),ISIperc_t600_t900.peak_mean = mean (ISIperc_t600_t900.peak, na.rm=TRUE), ISIperc_t600_t900.peak_median = median (ISIperc_t600_t900.peak, na.rm=TRUE),
                   ISIperc_t600_t900.peak_sd = sd(ISIperc_t600_t900.peak, na.rm=TRUE), ISIperc_t600_t900.peak_se = ISIperc_t600_t900.peak_sd/sqrt(ISIperc_t600_t900.peak_N), ISIperc_t600_t900.peak_cv = ISIperc_t600_t900.peak_sd/ISIperc_t600_t900.peak_mean, 
                   ISIperc_t600_t900.peak_CV2 = CV2 ((as.vector( ISIperc_t600_t900.peak ))[!is.na(as.vector( ISIperc_t600_t900.peak ))]),
                   
                   ISIperc_t650_t850.peak_N = sum (!is.na(ISIperc_t650_t850.peak)),ISIperc_t650_t850.peak_mean = mean (ISIperc_t650_t850.peak, na.rm=TRUE), ISIperc_t650_t850.peak_median = median (ISIperc_t650_t850.peak, na.rm=TRUE),
                   ISIperc_t650_t850.peak_sd = sd(ISIperc_t650_t850.peak, na.rm=TRUE), ISIperc_t650_t850.peak_se = ISIperc_t650_t850.peak_sd/sqrt(ISIperc_t650_t850.peak_N), ISIperc_t650_t850.peak_cv = ISIperc_t650_t850.peak_sd/ISIperc_t650_t850.peak_mean, 
                   ISIperc_t650_t850.peak_CV2 = CV2 ((as.vector(  ISIperc_t650_t850.peak ))[!is.na(as.vector(  ISIperc_t650_t850.peak ))]),
                   
                   ISIperc_t700_t800.peak_N = sum (!is.na(ISIperc_t700_t800.peak)),ISIperc_t700_t800.peak_mean = mean (ISIperc_t700_t800.peak, na.rm=TRUE), ISIperc_t700_t800.peak_median = median (ISIperc_t700_t800.peak, na.rm=TRUE),
                   ISIperc_t700_t800.peak_sd = sd(ISIperc_t700_t800.peak, na.rm=TRUE), ISIperc_t700_t800.peak_se = ISIperc_t700_t800.peak_sd/sqrt(ISIperc_t700_t800.peak_N), ISIperc_t700_t800.peak_cv = ISIperc_t700_t800.peak_sd/ISIperc_t700_t800.peak_mean, 
                   ISIperc_t700_t800.peak_CV2 = CV2 ((as.vector( ISIperc_t700_t800.peak ))[!is.na(as.vector( ISIperc_t700_t800.peak ))]),
                   
                   ISIperc_t725_t775.peak_N = sum (!is.na(ISIperc_t725_t775.peak)),ISIperc_t725_t775.peak_mean = mean (ISIperc_t725_t775.peak, na.rm=TRUE), ISIperc_t725_t775.peak_median = median (ISIperc_t725_t775.peak, na.rm=TRUE),
                   ISIperc_t725_t775.peak_sd = sd(ISIperc_t725_t775.peak, na.rm=TRUE), ISIperc_t725_t775.peak_se = ISIperc_t725_t775.peak_sd/sqrt(ISIperc_t725_t775.peak_N), ISIperc_t725_t775.peak_cv = ISIperc_t725_t775.peak_sd/ISIperc_t725_t775.peak_mean, 
                   ISIperc_t725_t775.peak_CV2 = CV2 ((as.vector(  ISIperc_t725_t775.peak ))[!is.na(as.vector(  ISIperc_t725_t775.peak ))]),
                   
                   ISIperc_t650_t1000.peak_N = sum (!is.na(ISIperc_t650_t1000.peak)),ISIperc_t650_t1000.peak_mean = mean (ISIperc_t650_t1000.peak, na.rm=TRUE), ISIperc_t650_t1000.peak_median = median (ISIperc_t650_t1000.peak, na.rm=TRUE),
                   ISIperc_t650_t1000.peak_sd = sd(ISIperc_t650_t1000.peak, na.rm=TRUE), ISIperc_t650_t1000.peak_se = ISIperc_t650_t1000.peak_sd/sqrt(ISIperc_t650_t1000.peak_N), ISIperc_t650_t1000.peak_cv = ISIperc_t650_t1000.peak_sd/ISIperc_t650_t1000.peak_mean, 
                   ISIperc_t650_t1000.peak_CV2 = CV2 ((as.vector(ISIperc_t650_t1000.peak  ))[!is.na(as.vector(ISIperc_t650_t1000.peak  ))]),
                   
                   
                   CRperc_t600_t900.hp_N = sum (!is.na(CRperc_t600_t900.hp)),CRperc_t600_t900.hp_mean = mean (CRperc_t600_t900.hp, na.rm=TRUE), CRperc_t600_t900.hp_median = median (CRperc_t600_t900.hp, na.rm=TRUE),
                   CRperc_t600_t900.hp_sd = sd(CRperc_t600_t900.hp, na.rm=TRUE), CRperc_t600_t900.hp_se = CRperc_t600_t900.hp_sd/sqrt(CRperc_t600_t900.hp_N), CRperc_t600_t900.hp_cv = CRperc_t600_t900.hp_sd/CRperc_t600_t900.hp_mean, 
                   CRperc_t600_t900.hp_CV2 = CV2 ((as.vector( CRperc_t600_t900.hp ))[!is.na(as.vector( CRperc_t600_t900.hp ))]),
                   
                   CRperc_t650_t850.hp_N = sum (!is.na(CRperc_t650_t850.hp)),CRperc_t650_t850.hp_mean = mean (CRperc_t650_t850.hp, na.rm=TRUE), CRperc_t650_t850.hp_median = median (CRperc_t650_t850.hp, na.rm=TRUE),
                   CRperc_t650_t850.hp_sd = sd(CRperc_t650_t850.hp, na.rm=TRUE), CRperc_t650_t850.hp_se = CRperc_t650_t850.hp_sd/sqrt(CRperc_t650_t850.hp_N), CRperc_t650_t850.hp_cv = CRperc_t650_t850.hp_sd/CRperc_t650_t850.hp_mean, 
                   CRperc_t650_t850.hp = CV2 ((as.vector( CRperc_t650_t850.hp ))[!is.na(as.vector( CRperc_t650_t850.hp ))]),
                   
                   CRperc_t700_t800.hp_N = sum (!is.na(CRperc_t700_t800.hp)),CRperc_t700_t800.hp_mean = mean (CRperc_t700_t800.hp, na.rm=TRUE), CRperc_t700_t800.hp_median = median (CRperc_t700_t800.hp, na.rm=TRUE),
                   CRperc_t700_t800.hp_sd = sd(CRperc_t700_t800.hp, na.rm=TRUE), CRperc_t700_t800.hp_se = CRperc_t700_t800.hp_sd/sqrt(CRperc_t700_t800.hp_N), CRperc_t700_t800.hp_cv = CRperc_t700_t800.hp_sd/CRperc_t700_t800.hp_mean, 
                   CRperc_t700_t800.hp_CV2 = CV2 ((as.vector( CRperc_t700_t800.hp ))[!is.na(as.vector( CRperc_t700_t800.hp ))]),
                   
                   CRperc_t725_t775.hp_N = sum (!is.na(CRperc_t725_t775.hp)),CRperc_t725_t775.hp_mean = mean (CRperc_t725_t775.hp, na.rm=TRUE), CRperc_t725_t775.hp_median = median (CRperc_t725_t775.hp, na.rm=TRUE),
                   CRperc_t725_t775.hp_sd = sd(CRperc_t725_t775.hp, na.rm=TRUE), CRperc_t725_t775.hp_se = CRperc_t725_t775.hp_sd/sqrt(CRperc_t725_t775.hp_N), CRperc_t725_t775.hp_cv = CRperc_t725_t775.hp_sd/CRperc_t725_t775.hp_mean, 
                   CRperc_t725_t775.hp_CV2 = CV2 ((as.vector(CRperc_t725_t775.hp  ))[!is.na(as.vector( CRperc_t725_t775.hp ))]),
                   
                   CRperc_t650_t1000.hp_N = sum (!is.na(CRperc_t650_t1000.hp)),CRperc_t650_t1000.hp_mean = mean (CRperc_t650_t1000.hp, na.rm=TRUE), CRperc_t650_t1000.hp_median = median (CRperc_t650_t1000.hp, na.rm=TRUE),
                   CRperc_t650_t1000.hp_sd = sd(CRperc_t650_t1000.hp, na.rm=TRUE), CRperc_t650_t1000.hp_se = CRperc_t650_t1000.hp_sd/sqrt(CRperc_t650_t1000.hp_N), CRperc_t650_t1000.hp_cv = CRperc_t650_t1000.hp_sd/CRperc_t650_t1000.hp_mean, 
                   CRperc_t650_t1000.hp_CV2 = CV2 ((as.vector(CRperc_t650_t1000.hp  ))[!is.na(as.vector( CRperc_t650_t1000.hp ))]),
                   
                   CRperc_t600_t900.peak_N = sum (!is.na(CRperc_t600_t900.peak)),CRperc_t600_t900.peak_mean = mean (CRperc_t600_t900.peak, na.rm=TRUE), CRperc_t600_t900.peak_median = median (CRperc_t600_t900.peak, na.rm=TRUE),
                   CRperc_t600_t900.peak_sd = sd(CRperc_t600_t900.peak, na.rm=TRUE), CRperc_t600_t900.peak_se = CRperc_t600_t900.peak_sd/sqrt(CRperc_t600_t900.peak_N), CRperc_t600_t900.peak_cv = CRperc_t600_t900.peak_sd/CRperc_t600_t900.peak_mean, 
                   CRperc_t600_t900.peak_CV2 = CV2 ((as.vector( CRperc_t600_t900.peak ))[!is.na(as.vector( CRperc_t600_t900.peak ))]),
                   
                   CRperc_t650_t850.peak_N = sum (!is.na(CRperc_t650_t850.peak)),CRperc_t650_t850.peak_mean = mean (CRperc_t650_t850.peak, na.rm=TRUE), CRperc_t650_t850.peak_median = median (CRperc_t650_t850.peak, na.rm=TRUE),
                   CRperc_t650_t850.peak_sd = sd(CRperc_t650_t850.peak, na.rm=TRUE), CRperc_t650_t850.peak_se = CRperc_t650_t850.peak_sd/sqrt(CRperc_t650_t850.peak_N), CRperc_t650_t850.peak_cv = CRperc_t650_t850.peak_sd/CRperc_t650_t850.peak_mean, 
                   CRperc_t650_t850.peak_CV2 = CV2 ((as.vector( CRperc_t650_t850.peak ))[!is.na(as.vector( CRperc_t650_t850.peak  ))]),
                   
                   CRperc_t700_t800.peak_N = sum (!is.na(CRperc_t700_t800.peak)),CRperc_t700_t800.peak_mean = mean (CRperc_t700_t800.peak, na.rm=TRUE), CRperc_t700_t800.peak_median = median (CRperc_t700_t800.peak, na.rm=TRUE),
                   CRperc_t700_t800.peak_sd = sd(CRperc_t700_t800.peak, na.rm=TRUE), CRperc_t700_t800.peak_se = CRperc_t700_t800.peak_sd/sqrt(CRperc_t700_t800.peak_N), CRperc_t700_t800.peak_cv = CRperc_t700_t800.peak_sd/CRperc_t700_t800.peak_mean, 
                   CRperc_t700_t800.peak_CV2 = CV2 ((as.vector(CRperc_t700_t800.peak  ))[!is.na(as.vector( CRperc_t700_t800.peak ))]),
                   
                   CRperc_t725_t775.peak_N = sum (!is.na(CRperc_t725_t775.peak)),CRperc_t725_t775.peak_mean = mean (CRperc_t725_t775.peak, na.rm=TRUE), CRperc_t725_t775.peak_median = median (CRperc_t725_t775.peak, na.rm=TRUE),
                   CRperc_t725_t775.peak_sd = sd(CRperc_t725_t775.peak, na.rm=TRUE), CRperc_t725_t775.peak_se = CRperc_t725_t775.peak_sd/sqrt(CRperc_t725_t775.peak_N), CRperc_t725_t775.peak_cv = CRperc_t725_t775.peak_sd/CRperc_t725_t775.peak_mean,
                   CRperc_t725_t775.peak_CV2 = CV2 ((as.vector( CRperc_t725_t775.peak ))[!is.na(as.vector( CRperc_t725_t775.peak ))]),
                   
                   CRperc_t650_t1000.peak_N = sum (!is.na(CRperc_t650_t1000.peak)),CRperc_t650_t1000.peak_mean = mean (CRperc_t650_t1000.peak, na.rm=TRUE), CRperc_t650_t1000.peak_median = median (CRperc_t650_t1000.peak, na.rm=TRUE),
                   CRperc_t650_t1000.peak_sd = sd(CRperc_t650_t1000.peak, na.rm=TRUE), CRperc_t650_t1000.peak_se = CRperc_t650_t1000.peak_sd/sqrt(CRperc_t650_t1000.peak_N), CRperc_t650_t1000.peak_cv = CRperc_t650_t1000.peak_sd/CRperc_t650_t1000.peak_mean, 
                   CRperc_t650_t1000.peak_CV2 = CV2 ((as.vector(CRperc_t650_t1000.peak  ))[!is.na(as.vector( CRperc_t650_t1000.peak ))])  )
}

### xxxxxxx_CV2 = CV2 ((as.vector(  ))[!is.na(as.vector(  ))])

# remove extreme outliers from CV and CV2 columns
for (u in 1:length(colnames(ebcAggr))) 
{
  col <- colnames(ebcAggr)[u]
  ebcAggr[[col]]
  if ((grepl("_CV2", col) == TRUE) | (grepl("_cv", col) == TRUE )) 
  {
    # boxplot (ebcAggr[[col]])
    # print(sort(ebcAggr[[col]]))
    # print(quantile(ebcAggr[[col]], na.rm=TRUE))
    # print(5*(quantile(ebcAggr[[col]], na.rm = TRUE)[4]))
    ebcAggr[[col]] <- ifelse(ebcAggr[[col]] <= (5*(quantile(ebcAggr[[col]], na.rm = TRUE)[4])), ebcAggr[[col]], NA )
    # boxplot (ebcAggr[[col]])
    # print(sort(ebcAggr[[col]]))
  } else {}
  
}

# define amplitude covariate
ebcAggr$amplitude_covariate <- as.numeric(ebcAggr$AmpattX_mean)

# recode session into numeric for models 1, 2, and 3.
ebcAggr$session_nr2 <- as.numeric(ebcAggr$session_nr)
ebcAggr$block_nr2 <- as.numeric(ebcAggr$block_nr)

# view ebcAggr
nrow(ebcAggr)
dim(ebcAggr)
write.csv(ebcAggr, file = "ebcAggr.csv")


# make aggregate file 5 which is aggregate of aggregate
rm(ebcAggr5)
{ebcAggr5 <- ebcAggr %>% group_by (genotype, session_nr) %>% summarise_if (is.numeric, funs( 
  mean (., na.rm=T),
  median (., na.rm=T),
  n = sum(!is.na(.)),
  #min = min (.,na.rm=T),
  #max = max (.,na.rm=T),
  sd = sd (.,na.rm=T),
  se = sd(., na.rm=T)/sqrt(sum(!is.na(.))), 
  cv = sd(., na.rm=T)/mean (., na.rm=T),
  # meanCI = ci (., na.rm=T)[1],
  lowCI = ci (., na.rm=T)[2],
  highCI = ci (., na.rm=T)[3],
  # sdCI = ci (., na.rm=T)[4]
  ci95 = ((ci (., na.rm=T)[3])-(ci (., na.rm=T)[2]))/2
))
}

ebcAggr5 <- ebcAggr5 [ , order(names(ebcAggr5))]
write.csv(ebcAggr5, file = "ebcAggr5.csv")

# make array with names of dataframes to be analyzed
dataframes <- c("ebcAll", "ebcAggr")



################################
##PERFORM LME ON ALL VARIABLES##
################################

for (y in 1:(length(dataframes)))
{
  #convert dataframe to name 'ebc'
  {rm(ebc)
  ebc <- get(dataframes [y])} #convert string to object
  #View(ebc)
  
  # make new folder in the wd and use it to store output
  {foldername <- (paste0("R-analysis folder ", dataframes [y]))
  dir.create(paste0 (wd, "/", foldername), showWarnings = FALSE)}
  
  # find column names in dataset
  {
    columnName <- colnames(ebc)
    print(columnName)
  }
  
  # set LME settings to optimal (more than 10 iterations etc)
  {
    ctrl <- lmeControl(opt='optim')
    # ctrl <- lmeControl(opt='optim', optimMethod = "SANN")
  }
  
  # analyze selected columns in dataframe using for loop. For alldata use column numbers 17:33, For aggr use column numbers 5:56
  {
  #for (i in 1:length(columnName))
  for (i in 1:length(columnName))
    tryCatch(
    
    # if (is.numeric(ebc[,columnName[i]]) == TRUE & (columnName[i] != "session_nr2")) # analyze all numeric columns
    # if (is.numeric(ebc[,columnName[i]]) == TRUE & (columnName[i] != "session_nr2") & ((grepl("time", columnName[i]) == TRUE) | (grepl("onset", columnName[i]) == TRUE))) # only timing data will be analyzed
    # if (is.numeric(ebc[,columnName[i]]) == TRUE & (columnName[i] != "session_nr2") & (grepl("mp", columnName[i]) == TRUE)) # only amplitude data will be analyzed
    if  (
        (columnName[i] == "CRperc5.hp_mean") |
        (columnName[i] == "CRperc10.hp_mean") |
        (columnName[i] == "CRperc5.peak_mean") |
        (columnName[i] == "CRperc10.peak_mean") |
        (columnName[i] == "CRperc_t650_t1000.hp_mean") |
        (columnName[i] == "CRperc_t700_t800.peak_mean") |

        (columnName[i] == "ISIpeakamp") |
        (columnName[i] == "ISIpeakamp_mean") |
        (columnName[i] == "ISIpeakamp_cv") |
        (columnName[i] == "ISIpeakamp_sd") |

        (columnName[i] == "log2ISIpeakamp") |
        (columnName[i] == "log2ISIpeakamp_mean") |
        (columnName[i] == "log2ISIpeakamp_cv") |
        (columnName[i] == "log2ISIpeakamp_sd") |

        (columnName[i] == "ISIpeaktime") |
        (columnName[i] == "ISIpeaktime_mean") |
        (columnName[i] == "ISIpeaktime_cv") |
        (columnName[i] == "ISIpeaktime_sd") |

        (columnName[i] == "ISIhpamp") |
        (columnName[i] == "ISIhpamp_mean") |
        (columnName[i] == "ISIhpamp_cv") |
        (columnName[i] == "ISIhpamp_sd") |

        (columnName[i] == "log2ISIhpamp") |
        (columnName[i] == "log2ISIhpamp_mean") |
        (columnName[i] == "log2ISIhpamp_cv") |
        (columnName[i] == "log2ISIhpamp_sd") |

        (columnName[i] == "ISIhptime") |
        (columnName[i] == "ISIhptime_mean") |
        (columnName[i] == "ISIhptime_cv") |
        (columnName[i] == "ISIhptime_sd") |

        (columnName[i] == "AmpattX") |
        (columnName[i] == "AmpattX_mean") |
        (columnName[i] == "AmpattX_cv") |
        (columnName[i] == "AmpattX_sd") |

        (columnName[i] == "log2AmpattX") |
        (columnName[i] == "log2AmpattX_mean") |
        (columnName[i] == "log2AmpattX_cv") |
        (columnName[i] == "log2AmpattX_sd") |

        (columnName[i] == "CRonset") |
        (columnName[i] == "CRonset_mean") |
        (columnName[i] == "CRonset_cv") |
        (columnName[i] == "CRonset_sd") |

        (columnName[i] == "CRpeakamp") |
        (columnName[i] == "CRpeakamp_mean") |
        (columnName[i] == "CRpeakamp_cv") |
        (columnName[i] == "CRpeakamp_sd") |

        (columnName[i] == "log2CRpeakamp") |
        (columnName[i] == "log2CRpeakamp_mean") |
        (columnName[i] == "log2CRpeakamp_cv") |
        (columnName[i] == "log2CRpeakamp_sd") |

        (columnName[i] == "CRpeaktime") |
        (columnName[i] == "CRpeaktime_mean") |
        (columnName[i] == "CRpeaktime_cv") |
        (columnName[i] == "CRpeaktime_sd") |

        (columnName[i] == "CRhptime") |
        (columnName[i] == "CRhptime_mean") |
        (columnName[i] == "CRhptime_cv") |
        (columnName[i] == "CRhptime_sd") |

        (columnName[i] == "CRhpamp") |
        (columnName[i] == "CRhpamp_mean") |
        (columnName[i] == "CRhpamp_cv") |
        (columnName[i] == "CRhpamp_sd")
        # 
        # (columnName[i] == "log2CRhpamp") |
        # (columnName[i] == "log2CRhpamp_mean") |
        # (columnName[i] == "log2CRhpamp_cv") |
        # (columnName[i] == "log2CRhpamp_sd") |
        # 
        # (columnName[i] == "CRpeakscore") |
        # (columnName[i] == "CRpeakscore_mean") |
        # (columnName[i] == "CRhpscore") |
        # (columnName[i] == "CRhpscore_mean")
        ) # analyze specific columns for acquisition sessions
    
    {
      {
          # show variable to be analyzed
          {
          variable <- columnName[i]
          print (variable)
          }
          
          # create new subfolder in R-analysis to save all output in textfiles
          {
          dir.create(paste0 (wd, "/", foldername, "/", variable, "_", i), showWarnings = FALSE)
          }
        
          # make plots
        if (unique (ebc$session_nr2) == 1) { 
          graph_pos <- "session_nr ~ groups"}
        else {
          graph_pos <- "groups.~"}
        
        
          {
            # make boxplot of data and export
            {boxplot <- ggplot(data=ebc, aes(x=session_nr, y=get(colnames(ebc)[i]), group=session_nr)) +
              geom_boxplot (aes(color=genotype)) +
              facet_grid(genotype~.) +
              scale_y_continuous(name=(colnames(ebc)[i])) + 
              theme_classic () +
              stat_boxplot(geom ='errorbar', width = 0.3, aes(color=genotype)) 
            print(boxplot)
            dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_boxplot", ".pdf")))
            while (dev.cur()>1) dev.off()
            if(!is.null(dev.list())) dev.off()}
            
            # make dotboxplot of data and export  
            {dotboxplot <- ggplot(data=ebc, aes(x=session_nr, y=get(colnames(ebc)[i]), group=session_nr)) +
                geom_boxplot(aes(color=genotype)) +
                geom_dotplot (binaxis='y', stackdir='center', dotsize=0.3) +
                facet_grid(genotype~.) +
                scale_y_continuous(name=(colnames(ebc)[i])) + 
                stat_boxplot(geom ='errorbar', width = 0.3, aes(color=genotype)) +
                theme_classic ()
              print(dotboxplot)
              dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_dotboxplot", ".pdf")))
              while (dev.cur()>1) dev.off()
              if(!is.null(dev.list())) dev.off()}
            
            # make lineplot only on aggregate date, histogram only on alldata
            if (grepl("All", dataframes [y]) == TRUE) 
              # make histogram and freqpoly plot
            {
              # determine binwidth
              {
                if (grepl("abs", variable) == TRUE) # all absolute values range from 1-10
                {(binw <- 1) & (xlim <- c(-2,10))
                } else if (grepl("perc", variable) == TRUE) # all percentage values range from 0-100
                {(binw <- 10) & (xlim <- c(-10,110))
                } else if (grepl("onset", variable) == TRUE) # all onset have normal timing range from 500 to 750
                {(binw <- 20) & (xlim <- c(400,1600))
                } else if ((grepl("amp", variable) == TRUE) & (grepl("abs", variable) == FALSE) & (grepl("norm", variable) == FALSE))
                {(binw <- 0.1) & (xlim <- c(-0.5,2))
                } else if ((grepl("amp", variable) == TRUE) & (grepl("abs", variable) == FALSE)) # real amplitude values range from -1 to 2
                {(binw <- 10) & (xlim <- c(-50,150))
                } else if ((grepl("time", variable) == TRUE) & (grepl("abs", variable) == FALSE))
                {(binw <- 20) & (xlim <- c(400,1600))
                } else if ((grepl("amp", variable) == TRUE) & (grepl("abs", variable) == FALSE))
                {(binw <- 10) & (xlim <- c(-50,150))
                } else {(binw <- NULL) & (xlim <- NULL)}
              }
              

              # make histogram and export
              {histogram1 <- ggplot(data=ebc, aes(get(colnames(ebc)[i]))) +
                  geom_histogram ((aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), binwidth=(binw)) +
                  facet_grid(.~ genotype) +
                  scale_x_continuous(name=(colnames(ebc)[i]), limits = (xlim)) +
                  scale_y_continuous(name="Frequency percentage/bin", limits = (c(0,100))) + 
                  theme_classic ()
                print(histogram1)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_histogram1", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              {histogram2 <- ggplot(data=ebc, aes(get(colnames(ebc)[i]))) +
                  geom_histogram ((aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), binwidth=(binw)) +
                  facet_grid(session_nr ~ genotype) +
                  scale_x_continuous(name=(colnames(ebc)[i]), limits = (xlim)) +
                  scale_y_continuous(name="Frequency percentage/bin", limits = (c(0,100))) + 
                  theme_classic ()
                print(histogram2)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_histogram2", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              
        
              # make freqpoly and export
              {freqpoly1 <- ggplot(data=ebc, aes(get(colnames(ebc)[i]))) +
                  geom_freqpoly ((aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), binwidth=(binw)) +
                  facet_grid(session_nr ~ genotype) +
                  scale_x_continuous(name = (colnames(ebc)[i]), limits = (xlim)) +
                  scale_y_continuous(name="Frequency percentage/bin", limits = (c(0,100))) +
                  theme_classic ()
                print(freqpoly1)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_freqpoly1", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              {freqpoly2 <- ggplot(data=ebc, aes(get(colnames(ebc)[i]))) +
                  geom_freqpoly ((aes(y = (..count..)/tapply(..count..,..PANEL..,sum)[..PANEL..]*100)), binwidth=(binw)) +
                  facet_grid(.~genotype) +
                  scale_x_continuous(name = (colnames(ebc)[i]), limits = (xlim)) +
                  scale_y_continuous(name="Frequency percentage/bin", limits = (c(0,100))) +
                  theme_classic ()
                print(freqpoly2)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_freqpoly2", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make density with bars
              {density1 <- ggplot(data=ebc, aes(get(colnames(ebc)[i]), fill = genotype, color = genotype)) +
                  facet_grid (.~genotype) +
                  theme_classic () +
                  geom_density () +
                  geom_histogram (aes (y = ..density..), color = "black", fill= "white", binwidth = (binw)) +
                  scale_x_continuous(name = (colnames(ebc)[i]), limits = (xlim)) +
                  scale_y_continuous(name="Density/bin")
                print(density1)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_density1", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              {density2 <- ggplot(data=ebc, aes(get(colnames(ebc)[i]), fill = genotype, color = genotype)) +
                  facet_grid (session_nr ~ genotype) +
                  theme_classic () +
                  geom_density () +
                  geom_histogram (aes (y = ..density..), color = "black", fill= "white", binwidth = (binw)) +
                  scale_x_continuous(name = (colnames(ebc)[i]), limits = (xlim)) +
                  scale_y_continuous(name="Density/bin")
                print(density2)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_density2", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make csv file with mean sd se ci and export
              {
                rm(ebcAggr3)
                ebcAggr3 <- (data.frame(summarySE(ebc, measurevar = variable, groupvar = c("genotype", "session_nr2", "block_nr2"), na.rm = TRUE)))
                outLME <- capture.output(ebcAggr3)
                cat (paste0("Title:", variable, "_Aggr3"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_Aggr3", ".txt")), sep = "\n", append = FALSE)}
              
              # make trial-by-trial plot
              {trialbytrial <- ggplot(data = ebcAggr3, aes(x=block_nr2, y=(get(colnames(ebc)[i])), group=genotype)) +
                  geom_point (aes(color = genotype)) +
                  geom_smooth(method=loess, se = FALSE) +
                  facet_grid(genotype ~ session_nr2) +
                  scale_x_continuous(limits = c(0,20)) +
                  scale_y_continuous(name = (colnames(ebc)[i]), limits = (xlim)) +
                  theme_classic()
                print(trialbytrial)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_trialbytrial", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()
                
              }
              
              
            }
            
            else if (grepl("Aggr", dataframes [y]) == TRUE) 
              # make mean and median lineplot
            {
              # make mean lineplot of data and export
              {
                lineplot_mean <- ggplot(data=ebc, aes(x=session_nr, y=get(colnames(ebc)[i]), group=mouse_id)) +
                  geom_line(aes(color=genotype)) +
                  geom_point(aes(color=genotype)) +
                  facet_grid(.~genotype) +
                  ggtitle (paste0 ("lineplot_mean_", colnames(ebc)[i])) +
                  stat_summary(fun = "mean", geom = "line", color = "black", size = 1.2, group = "genotype") +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme_classic() +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_mean)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_mean", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make median lineplot of data and export
              {
                lineplot_median <- ggplot(data=ebc, aes(x=session_nr, y=get(colnames(ebc)[i]), group=mouse_id)) +
                  geom_line(aes(color=genotype)) +
                  geom_point(aes(color=genotype)) +
                  facet_grid(.~genotype) +
                  ggtitle (paste0 ("lineplot_median_", colnames(ebc)[i])) +
                  stat_summary(fun = "median", geom = "line", color = "black", size = 1.2, group = "genotype") +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme_classic() +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_median)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_median", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()} 
              
              # make csv file with mean sd se ci and export
              {
                rm(ebcAggr2)
                ebcAggr2 <- (data.frame(summarySE(ebcAggr, measurevar = variable, groupvar = c("genotype", "session_nr"), na.rm = TRUE)))
                outLME <- capture.output(ebcAggr2)
                # cat (paste0("Title:", variable, "_Aggr2"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_Aggr2", ".txt")), sep = "\n", append = FALSE)}
                write.csv(ebcAggr2, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_Aggr2.csv")))}
                
              # make lineplot_mean_se1 and export
              {lineplot_mean_se1 <- ggplot(data=ebcAggr2, aes(x=session_nr, y=get(colnames(ebc)[i]), group=genotype, color = genotype)) +
                  geom_errorbar(data = ebcAggr2, aes(ymin=get(colnames(ebc)[i]) - se, ymax = get(colnames(ebc)[i]) + se), width = .5, size = 1) +  
                  geom_line(size = 1) +
                  geom_point(size = 3) +
                  ggtitle (paste0 ("lineplot_mean_se1_", colnames(ebc)[i])) +
                  theme_classic() +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_mean_se1)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_se1_mean", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make lineplot_mean_se2 and export
              {lineplot_mean_se2 <- ggplot(data=ebcAggr2, aes(x=session_nr, y=get(colnames(ebc)[i]), group=genotype, color = genotype)) +
                  geom_errorbar(data = ebcAggr2, aes(ymin=get(colnames(ebc)[i]) - se, ymax = get(colnames(ebc)[i]) + se), width = .5, size = 1) +  
                  geom_line(size = 1) +
                  geom_point(size = 2) +
                  ggtitle (paste0 ("lineplot_mean_se2_", colnames(ebc)[i])) +
                  facet_grid(.~ genotype) +
                  theme_classic() +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_mean_se2)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_se2_mean", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make lineplot_mean_sd1 and export
              {lineplot_mean_sd1 <- ggplot(data=ebcAggr2, aes(x=session_nr, y=get(colnames(ebc)[i]), group=genotype, color = genotype)) +
                  geom_errorbar(data = ebcAggr2, aes(ymin=get(colnames(ebc)[i]) - sd, ymax = get(colnames(ebc)[i]) + sd), width = .5, size = 1) +  
                  geom_line(size = 1) +
                  geom_point(size = 3) +
                  ggtitle (paste0 ("lineplot_mean_sd1_", colnames(ebc)[i])) +
                  theme_classic() +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_mean_sd1)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_sd1_mean", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make lineplot_mean_sd2 and export
              {lineplot_mean_sd2 <- ggplot(data=ebcAggr2, aes(x=session_nr, y=get(colnames(ebc)[i]), group=genotype, color = genotype)) +
                  geom_errorbar(data = ebcAggr2, aes(ymin=get(colnames(ebc)[i]) - sd, ymax = get(colnames(ebc)[i]) + sd), width = .5, size = 1) +  
                  geom_line(size = 1) +
                  geom_point(size = 2) +
                  ggtitle (paste0 ("lineplot_mean_sd2_", colnames(ebc)[i])) +
                  facet_grid(.~ genotype) +
                  theme_classic() +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_mean_sd2)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_sd2_mean", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make lineplot_mean_ci1 and export
              {lineplot_mean_ci1 <- ggplot(data=ebcAggr2, aes(x=session_nr, y=get(colnames(ebc)[i]), group=genotype, color = genotype)) +
                  geom_errorbar(data = ebcAggr2, aes(ymin=get(colnames(ebc)[i]) - ci, ymax = get(colnames(ebc)[i]) + ci), width = .5, size = 1) +  
                  geom_line(size = 1) +
                  geom_point(size = 3) +
                  theme_classic() +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_mean_ci1)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_ci1_mean", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
              
              # make lineplot_mean_ci2 and export
              {lineplot_mean_ci2 <- ggplot(data=ebcAggr2, aes(x=session_nr, y=get(colnames(ebc)[i]), group=genotype, color = genotype)) +
                  geom_errorbar(data = ebcAggr2, aes(ymin=get(colnames(ebc)[i]) - ci, ymax = get(colnames(ebc)[i]) + ci), width = .5, size = 1) +  
                  geom_line(size = 1) +
                  geom_point(size = 2) +
                  facet_grid(.~ genotype) +
                  theme_classic() +
                  scale_y_continuous(colnames(ebc)[i]) +
                  theme(axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                        axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                        axis.line = element_line(colour = 'black', size = 1), 
                        axis.ticks = element_line(colour = "black", size = 1))
                print(lineplot_mean_ci2)
                dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_lineplot_ci2_mean", ".pdf")))
                while (dev.cur()>1) dev.off()
                if(!is.null(dev.list())) dev.off()}
             
            }
            
            else {}
        }  
          
          # define model names
          {
          LME_variable0 <- paste0 ("LME_", variable, "0")
          LME_variable0c <- paste0 ("LME_", variable, "0_corr")
          LME_variable1 <- paste0 ("LME_", variable, "1")
          LME_variable1c <- paste0 ("LME_", variable, "1_corr")
          LME_variable2 <- paste0 ("LME_", variable, "2")
          LME_variable2c <- paste0 ("LME_", variable, "2_corr")
          LME_variable3 <- paste0 ("LME_", variable, "3")
          LME_variable3c <- paste0 ("LME_", variable, "3_corr")
          }
          # tryCatch(    
              # if there are more than 1 sessions, test model 0, 1, and 2.
              if (length (unique (ebc$session_nr2)) > 1) 
                {
                  # model 0: compound symmetry, RM anova equivalent. No correlation between sessions. Session_nr as factor.
                  LME_variable0 <-  lme (formula (paste ((colnames(ebc)[i]), "~ genotype + session_nr + session_nr:genotype")),
                                    control = ctrl, data = ebc, correlation = NULL, random = ~ 1 | mouse_id, method = "REML", na.action=na.exclude)
                  
                    # show results of model 0 and export
                    {
                        sum <- (summary (LME_variable0))
                        an <- (anova (LME_variable0))
                        vc <- (VarCorr (LME_variable0))
                        postHocGroupNoAdj <- emmeans (LME_variable0, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                        postHocGroupPerSesNoAdj <- (emmeans (LME_variable0, list (pairwise ~ genotype | session_nr), adjust = "no")) # posthoc test: group per session effect with p value
                        postHocGroupHolmAdj <- emmeans (LME_variable0, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                        postHocGroupPerSesHolmAdj <- (emmeans (LME_variable0, list (pairwise ~ genotype | session_nr), adjust = "holm")) # posthoc test: group per session effect with p value
                        postHocGroupFDRAdj <- emmeans (LME_variable0, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                        postHocGroupPerSesFDRAdj <- (emmeans (LME_variable0, list (pairwise ~ genotype | session_nr), adjust = "FDR")) # posthoc test: group per session effect with p value
                        
                      #   combs <- combn(4, 2)
                      #   p_values <- numeric(ncol(combs))
                      #   for (j in 1:ncol(combs)) {
                      #   if (combs[1, j] == 1) {
                      #     nams1 <- (names(fixef(LME_variable0)))
                      #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                      #     nam_group <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                      #     p_values[j] <- anova(LME_variable0, L = L)$'p-value'      
                      #   } else {
                      #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                      #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                      #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                      #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                      #     p_values[j] <- anova(LME_variable0, L = L)$'p-value'      
                      #   }
                      # }
                      #   postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                      #   postHocgenotypees
                      #   p_values
                  
                      
                      # export all output of LME_variable0 to textfile
                        outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupPerSesNoAdj, postHocGroupHolmAdj, postHocGroupPerSesHolmAdj, postHocGroupFDRAdj, postHocGroupPerSesFDRAdj)
                        cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model0", ".txt")), sep = "\n", append = FALSE)
                    }
                    
                    # plot model 0 and export
                    # { 
                    # # effect plot
                    #   if (length (unique (ebc$session_nr2)) >= 1) { # effect plot only shown if there is more than 1 session.
                    #     effectPlotData <- function (object, newdata, orig_data) {
                    #     form <- formula(object)
                    #     namesVars <- all.vars(form)
                    #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #     Terms <- delete.response(terms(form))
                    #     mfX <- model.frame(Terms, data = orig_data)
                    #     Terms_new <- attr(mfX, "terms")
                    #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #     X <- model.matrix(Terms_new, mfX_new)
                    #     pred <- c(X %*% betas)
                    #     ses <- sqrt(diag(X %*% V %*% t(X)))
                    #     newdata$pred <- pred
                    #     newdata$low <- pred - 1.96 * ses
                    #     newdata$upp <- pred + 1.96 * ses
                    #     newdata
                    #   }
                    #   newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr)))
                    #   plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                    #                   data = effectPlotData(LME_variable0, newDF, ebc),
                    #                   #lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                    #                   lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                    #                   xlab = "Session",
                    #                   ylab = (paste0 (variable)))
                    #   print(plot1)
                    # 
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model0", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    #   } else {
                    #   }
                    #   # Residual plot
                    #   rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                    #     form <- formula(object)
                    #     namesVars <- all.vars(form)
                    #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #     Terms <- delete.response(terms(form))
                    #     mfX <- model.frame(Terms, data = orig_data)
                    #     Terms_new <- attr(mfX, "terms")
                    #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #     X <- model.matrix(Terms_new, mfX_new)
                    #     ind <- combn(nrow(X), 2)
                    #     k <- ncol(ind)
                    #     out <- matrix(0, k, 5)
                    #     for (i in seq_len(k)) {
                    #       XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                    #       est <- drop(XX %*% betas)
                    #       se <- sqrt(diag(XX %*% V %*% t(XX)))
                    #       out[i, 1] <- est
                    #       out[i, 2] <- se
                    #       out[i, 3] <- est - 1.96 * se
                    #       out[i, 4] <- est + 1.96 * se
                    #       out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                    #     }
                    #     if (k > 2 && adjust.p) {
                    #       out[, 5] <- p.adjust(out[, 5], ...)
                    #     }
                    #     colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                    #     rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                    #     out
                    #   }
                    # 
                    #   v_pVal <- vector(mode="numeric", length=0)
                    #   v_outRowDiff <- vector(mode="numeric", length=0)
                    # 
                    #   for (k in 1:max(ebc$session_nr2))
                    #   {
                    #     # test the effect at specific sessions
                    #     nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype)))
                    #     outRowDiff <- capture.output(nDF, rowDiff(LME_variable0, nDF, ebc))
                    #     v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                    #     pVal <- rowDiff(LME_variable0, nDF, ebc) [5]
                    #     v_pVal <- append (v_pVal, pVal)
                    #   }
                    # 
                    #   v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                    # 
                    #   cat (paste0("Title:", variable, "_session_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_session_diff", "_model0", ".txt")), sep = "\n", append = FALSE)
                    # 
                    #   # plot residuals
                    #   plot2 <- plot(LME_variable0)
                    #   print(plot2)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model0", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # 
                    #   plot3 <- qqnorm(LME_variable0)
                    #   print(plot3)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model0", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # }
                   
                  
                  # if 'time' is in columnheader, test same model 0 but now use AmpattX as covariate.  
                  if ((grepl("time", variable) == TRUE) | (grepl("onset", variable) == TRUE)) {
                  # Model 0c: idem as model 0 but now with AmpattX as fixed factor added
                  LME_variable0c <- lme (formula (paste ((colnames(ebc)[i]), "~ genotype + session_nr + session_nr:genotype + amplitude_covariate")),
                                    control = ctrl, data = ebc, correlation = NULL, random = ~ 1 | mouse_id, method = "REML", na.action=na.exclude) 
                     
                    # show results of model 0c and export
                    {
                        sum <- (summary (LME_variable0c))
                        an <- (anova (LME_variable0c))
                        vc <- (VarCorr (LME_variable0c))
                        postHocGroupNoAdj <- emmeans (LME_variable0c, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                        postHocGroupPerSesNoAdj <- (emmeans (LME_variable0c, list (pairwise ~ genotype | session_nr), adjust = "no")) # posthoc test: group per session effect with p value
                        postHocGroupHolmAdj <- emmeans (LME_variable0c, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                        postHocGroupPerSesHolmAdj <- (emmeans (LME_variable0c, list (pairwise ~ genotype | session_nr), adjust = "holm")) # posthoc test: group per session effect with p value
                        postHocGroupFDRAdj <- emmeans (LME_variable0c, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                        postHocGroupPerSesFDRAdj <- (emmeans (LME_variable0c, list (pairwise ~ genotype | session_nr), adjust = "FDR")) # posthoc test: group per session effect with p value
                        
                        # combs <- combn(4, 2)
                        # p_values <- numeric(ncol(combs))
                        # for (j in 1:ncol(combs)) {
                        #   if (combs[1, j] == 1) {
                        #     nams1 <- (names(fixef(LME_variable0c)))
                        #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                        #     nam_group <- paste0("genotype", combs[2, j], ":")
                        #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                        #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                        #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                        #     p_values[j] <- anova(LME_variable0c, L = L)$'p-value'      
                        #   } else {
                        #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                        #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                        #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                        #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                        #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                        #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                        #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                        #     p_values[j] <- anova(LME_variable0c, L = L)$'p-value'      
                        #   }
                        # }
                        # 
                        # 
                        # postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                        #                          
                                                 
                        # export all output of LME_variable0c to textfile
                        outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupPerSesNoAdj, postHocGroupHolmAdj, postHocGroupPerSesHolmAdj, postHocGroupFDRAdj, postHocGroupPerSesFDRAdj)
                        cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model0c", ".txt")), sep = "\n", append = FALSE)
                      }
                      
                    # plot model 0c and export
                    # { 
                    # # effect plot
                    #     if (length (unique (ebc$session_nr2)) >= 1) { # effect plot only shown if there is more than 1 session.
                    #     effectPlotData <- function (object, newdata, orig_data) {
                    #       form <- formula(object)
                    #       namesVars <- all.vars(form)
                    #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #       Terms <- delete.response(terms(form))
                    #       mfX <- model.frame(Terms, data = orig_data)
                    #       Terms_new <- attr(mfX, "terms")
                    #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #       X <- model.matrix(Terms_new, mfX_new)
                    #       pred <- c(X %*% betas)
                    #       ses <- sqrt(diag(X %*% V %*% t(X)))
                    #       newdata$pred <- pred
                    #       newdata$low <- pred - 1.96 * ses
                    #       newdata$upp <- pred + 1.96 * ses
                    #       newdata
                    #     }
                    #     newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr), amplitude_covariate = median (amplitude_covariate) ))
                    #     plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                    #                     data = effectPlotData(LME_variable0c, newDF, ebc),
                    #                     lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                    #                     xlab = "Session",
                    #                     ylab = (paste0 (variable)))
                    #     print(plot1)
                    #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model0c", ".pdf")))
                    #     while (dev.cur()>1) dev.off()
                    #     if(!is.null(dev.list())) dev.off()
                    #     } else {
                    # 
                    #     }
                    # 
                    #     # Residual plot
                    #     rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                    #       form <- formula(object)
                    #       namesVars <- all.vars(form)
                    #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #       Terms <- delete.response(terms(form))
                    #       mfX <- model.frame(Terms, data = orig_data)
                    #       Terms_new <- attr(mfX, "terms")
                    #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #       X <- model.matrix(Terms_new, mfX_new)
                    #       ind <- combn(nrow(X), 2)
                    #       k <- ncol(ind)
                    #       out <- matrix(0, k, 5)
                    #       for (i in seq_len(k)) {
                    #         XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                    #         est <- drop(XX %*% betas)
                    #         se <- sqrt(diag(XX %*% V %*% t(XX)))
                    #         out[i, 1] <- est
                    #         out[i, 2] <- se
                    #         out[i, 3] <- est - 1.96 * se
                    #         out[i, 4] <- est + 1.96 * se
                    #         out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                    #       }
                    #       if (k > 2 && adjust.p) {
                    #         out[, 5] <- p.adjust(out[, 5], ...)
                    #       }
                    #       colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                    #       rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                    #       out
                    #     }
                    # 
                    #     v_pVal <- vector(mode="numeric", length=0)
                    #     v_outRowDiff <- vector(mode="numeric", length=0)
                    # 
                    #     for (k in 1:max(ebc$session_nr2))
                    #     {
                    #       # test the effect at specific sessions
                    #       nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype), amplitude_covariate = median (amplitude_covariate) ))
                    #       outRowDiff <- capture.output(nDF, rowDiff(LME_variable0c, nDF, ebc))
                    #       v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                    #       pVal <- rowDiff(LME_variable0c, nDF, ebc) [5]
                    #       v_pVal <- append (v_pVal, pVal)
                    #     }
                    # 
                    #     v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                    # 
                    #     cat (paste0("Title:", variable, "_session_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_session_diff", "_model0c", ".txt")), sep = "\n", append = FALSE)
                    # 
                    #     # plot residuals
                    #     plot2 <- plot(LME_variable0c)
                    #     print(plot2)
                    #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model0c", ".pdf")))
                    #     while (dev.cur()>1) dev.off()
                    #     if(!is.null(dev.list())) dev.off()
                    # 
                    #     plot3 <- qqnorm(LME_variable0c)
                    #     print(plot3)
                    #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model0c", ".pdf")))
                    #     while (dev.cur()>1) dev.off()
                    #     if(!is.null(dev.list())) dev.off()
                    # }
                  
                  } 
                  else {
                    
                  }
                  
                  
                  # Model 1: correlations between amplitude values change over time in a linear fashion. For instance, amplitude on session 1
                  # is stronger correlated to amplitude at session 2 than to amplitude at session 10.
                  LME_variable1 <-  lme (formula (paste ((colnames(ebc)[i]), "~ genotype + session_nr + session_nr:genotype")), 
                                    control = ctrl, data = ebc, random = ~ session_nr2 | mouse_id, method = "REML", na.action=na.exclude) 
                  
                    # show results of model 1 and export
                    {
                      sum <- (summary (LME_variable1))
                      an <- (anova (LME_variable1))
                      vc <- (VarCorr (LME_variable1))
                      postHocGroupNoAdj <- emmeans (LME_variable1, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                      postHocGroupPerSesNoAdj <- (emmeans (LME_variable1, list (pairwise ~ genotype | session_nr), adjust = "no")) # posthoc test: group per session effect with p value
                      postHocGroupHolmAdj <- emmeans (LME_variable1, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                      postHocGroupPerSesHolmAdj <- (emmeans (LME_variable1, list (pairwise ~ genotype | session_nr), adjust = "holm")) # posthoc test: group per session effect with p value
                      postHocGroupFDRAdj <- emmeans (LME_variable1, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                      postHocGroupPerSesFDRAdj <- (emmeans (LME_variable1, list (pairwise ~ genotype | session_nr), adjust = "FDR")) # posthoc test: group per session effect with p value
                     
                      # combs <- combn(4, 2)
                      # p_values <- numeric(ncol(combs))
                      # for (j in 1:ncol(combs)) {
                      #   if (combs[1, j] == 1) {
                      #     nams1 <- (names(fixef(LME_variable1)))
                      #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                      #     nam_group <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                      #     p_values[j] <- anova(LME_variable1, L = L)$'p-value'      
                      #   } else {
                      #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                      #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                      #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                      #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                      #     p_values[j] <- anova(LME_variable1, L = L)$'p-value'      
                      #   }
                      # }
                      # postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                      # 
                      # export all output of LME_variable1 to textfile
                      outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupPerSesNoAdj, postHocGroupHolmAdj, postHocGroupPerSesHolmAdj, postHocGroupFDRAdj, postHocGroupPerSesFDRAdj)
                      cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model1", ".txt")), sep = "\n", append = FALSE)
                    }
                    
                    # plot model 1 and export
                    # { 
                    # # effect plot
                    #   if (length (unique (ebc$session_nr2)) >= 1) { # effect plot only shown if there is more than 1 session.
                    #     effectPlotData <- function (object, newdata, orig_data) {
                    #       form <- formula(object)
                    #       namesVars <- all.vars(form)
                    #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #       Terms <- delete.response(terms(form))
                    #       mfX <- model.frame(Terms, data = orig_data)
                    #       Terms_new <- attr(mfX, "terms")
                    #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #       X <- model.matrix(Terms_new, mfX_new)
                    #       pred <- c(X %*% betas)
                    #       ses <- sqrt(diag(X %*% V %*% t(X)))
                    #       newdata$pred <- pred
                    #       newdata$low <- pred - 1.96 * ses
                    #       newdata$upp <- pred + 1.96 * ses
                    #       newdata
                    #     }
                    #     newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr)))
                    #     plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                    #                     data = effectPlotData(LME_variable1, newDF, ebc),
                    #                     lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                    #                     xlab = "Session",
                    #                     ylab = (paste0 (variable)))
                    #     print(plot1)
                    #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model1", ".pdf")))
                    #     while (dev.cur()>1) dev.off()
                    #     if(!is.null(dev.list())) dev.off()
                    #   } else {
                    #   }
                    #   # Residual plot
                    #   rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                    #     form <- formula(object)
                    #     namesVars <- all.vars(form)
                    #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #     Terms <- delete.response(terms(form))
                    #     mfX <- model.frame(Terms, data = orig_data)
                    #     Terms_new <- attr(mfX, "terms")
                    #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #     X <- model.matrix(Terms_new, mfX_new)
                    #     ind <- combn(nrow(X), 2)
                    #     k <- ncol(ind)
                    #     out <- matrix(0, k, 5)
                    #     for (i in seq_len(k)) {
                    #       XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                    #       est <- drop(XX %*% betas)
                    #       se <- sqrt(diag(XX %*% V %*% t(XX)))
                    #       out[i, 1] <- est
                    #       out[i, 2] <- se
                    #       out[i, 3] <- est - 1.96 * se
                    #       out[i, 4] <- est + 1.96 * se
                    #       out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                    #     }
                    #     if (k > 2 && adjust.p) {
                    #       out[, 5] <- p.adjust(out[, 5], ...)
                    #     }
                    #     colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                    #     rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                    #     out
                    #   }
                    # 
                    #   v_pVal <- vector(mode="numeric", length=0)
                    #   v_outRowDiff <- vector(mode="numeric", length=0)
                    # 
                    #   for (k in 1:max(ebc$session_nr2))
                    #   {
                    #     # test the effect at specific sessions
                    #     nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype)))
                    #     outRowDiff <- capture.output(nDF, rowDiff(LME_variable1, nDF, ebc))
                    #     v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                    #     pVal <- rowDiff(LME_variable1, nDF, ebc) [5]
                    #     v_pVal <- append (v_pVal, pVal)
                    #   }
                    # 
                    #   v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                    # 
                    #   cat (paste0("Title:", variable, "_session_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_session_diff", "_model1", ".txt")), sep = "\n", append = FALSE)
                    # 
                    #   # plot residuals
                    #   plot2 <- plot(LME_variable1)
                    #   print(plot2)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model1", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # 
                    #   plot3 <- qqnorm(LME_variable1)
                    #   print(plot3)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model1", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # }
                  
                  if ((grepl("time", variable) == TRUE) | (grepl("onset", variable) == TRUE)) {
                  # Model 1c idem to model 0 but now with amplitude added as fixed factor
                  LME_variable1c <-  lme (formula (paste ((colnames(ebc)[i]), "~ genotype + session_nr + session_nr:genotype + amplitude_covariate")),
                                    control = ctrl, data = ebc, random = ~ session_nr2 | mouse_id, method = "REML", na.action=na.exclude)
                
                    # show results of model 1c and export
                    {
                      sum <- (summary (LME_variable1c))
                      an <- (anova (LME_variable1c))
                      vc <- (VarCorr (LME_variable1c))
                      postHocGroupNoAdj <- emmeans (LME_variable1c, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                      postHocGroupPerSesNoAdj <- (emmeans (LME_variable1c, list (pairwise ~ genotype | session_nr), adjust = "no")) # posthoc test: group per session effect with p value
                      postHocGroupHolmAdj <- emmeans (LME_variable1c, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                      postHocGroupPerSesHolmAdj <- (emmeans (LME_variable1c, list (pairwise ~ genotype | session_nr), adjust = "holm")) # posthoc test: group per session effect with p value
                      postHocGroupFDRAdj <- emmeans (LME_variable1c, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                      postHocGroupPerSesFDRAdj <- (emmeans (LME_variable1c, list (pairwise ~ genotype | session_nr), adjust = "FDR")) # posthoc test: group per session effect with p value
                      
                      # combs <- combn(4, 2)
                      # p_values <- numeric(ncol(combs))
                      # for (j in 1:ncol(combs)) {
                      #   if (combs[1, j] == 1) {
                      #     nams1 <- (names(fixef(LME_variable1c)))
                      #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                      #     nam_group <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                      #     p_values[j] <- anova(LME_variable1c, L = L)$'p-value'      
                      #   } else {
                      #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                      #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                      #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                      #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                      #     p_values[j] <- anova(LME_variable1c, L = L)$'p-value'      
                      #   }
                      # }
                      # postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                      # 
                      # export all output of LME_variable1c to textfile
                      outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupPerSesNoAdj, postHocGroupHolmAdj, postHocGroupPerSesHolmAdj, postHocGroupFDRAdj, postHocGroupPerSesFDRAdj)
                      cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model1c", ".txt")), sep = "\n", append = FALSE)
                    }
                    
                    # plot model 1c and export
                    # { 
                    # # effect plot
                    #   if (length (unique (ebc$session_nr2)) >= 1) { # effect plot only shown if there is more than 1 session.
                    #     effectPlotData <- function (object, newdata, orig_data) {
                    #       form <- formula(object)
                    #       namesVars <- all.vars(form)
                    #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #       Terms <- delete.response(terms(form))
                    #       mfX <- model.frame(Terms, data = orig_data)
                    #       Terms_new <- attr(mfX, "terms")
                    #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #       X <- model.matrix(Terms_new, mfX_new)
                    #       pred <- c(X %*% betas)
                    #       ses <- sqrt(diag(X %*% V %*% t(X)))
                    #       newdata$pred <- pred
                    #       newdata$low <- pred - 1.96 * ses
                    #       newdata$upp <- pred + 1.96 * ses
                    #       newdata
                    #     }
                    #     newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr), amplitude_covariate = median (amplitude_covariate) ))
                    #     plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                    #                     data = effectPlotData(LME_variable1c, newDF, ebc),
                    #                     lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                    #                     xlab = "Session",
                    #                     ylab = (paste0 (variable)))
                    #     print(plot1)
                    #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model1c", ".pdf")))
                    #     while (dev.cur()>1) dev.off()
                    #     if(!is.null(dev.list())) dev.off()
                    #   } else {
                    # 
                    #   }
                    # 
                    #   # Residual plot
                    #   rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                    #     form <- formula(object)
                    #     namesVars <- all.vars(form)
                    #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #     Terms <- delete.response(terms(form))
                    #     mfX <- model.frame(Terms, data = orig_data)
                    #     Terms_new <- attr(mfX, "terms")
                    #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #     X <- model.matrix(Terms_new, mfX_new)
                    #     ind <- combn(nrow(X), 2)
                    #     k <- ncol(ind)
                    #     out <- matrix(0, k, 5)
                    #     for (i in seq_len(k)) {
                    #       XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                    #       est <- drop(XX %*% betas)
                    #       se <- sqrt(diag(XX %*% V %*% t(XX)))
                    #       out[i, 1] <- est
                    #       out[i, 2] <- se
                    #       out[i, 3] <- est - 1.96 * se
                    #       out[i, 4] <- est + 1.96 * se
                    #       out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                    #     }
                    #     if (k > 2 && adjust.p) {
                    #       out[, 5] <- p.adjust(out[, 5], ...)
                    #     }
                    #     colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                    #     rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                    #     out
                    #   }
                    # 
                    #   v_pVal <- vector(mode="numeric", length=0)
                    #   v_outRowDiff <- vector(mode="numeric", length=0)
                    # 
                    #   for (k in 1:max(ebc$session_nr2))
                    #   {
                    #     # test the effect at specific sessions
                    #     nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype), amplitude_covariate = median (amplitude_covariate) ))
                    #     outRowDiff <- capture.output(nDF, rowDiff(LME_variable1c, nDF, ebc))
                    #     v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                    #     pVal <- rowDiff(LME_variable1c, nDF, ebc) [5]
                    #     v_pVal <- append (v_pVal, pVal)
                    #   }
                    # 
                    #   v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                    # 
                    #   cat (paste0("Title:", variable, "_session_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_session_diff", "_model1c", ".txt")), sep = "\n", append = FALSE)
                    # 
                    #   # plot residuals
                    #   plot2 <- plot(LME_variable1c)
                    #   print(plot2)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model1c", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # 
                    #   plot3 <- qqnorm(LME_variable1c)
                    #   print(plot3)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model1c", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # } 
                  
                  } else {}
                        
                  # Model 2: correlations between amplitude values change over time in a sqr fashion. For instance, amplitude on session 1
                  # is stronger correlated to amplitude at session 2 than to amplitude at session 10.
                  LME_variable2 <-  lme (formula (paste ((colnames(ebc)[i]), "~ genotype + session_nr + session_nr:genotype")),
                                    control = ctrl, data = ebc, random = list(mouse_id = pdDiag(form = ~ poly(session_nr2, 2))), method = "REML", na.action=na.exclude) 
                  
                    # show results of model 2 and export
                    {
                      sum <- (summary (LME_variable2))
                      an <- (anova (LME_variable2))
                      vc <- (VarCorr (LME_variable2))
                      postHocGroupNoAdj <- emmeans (LME_variable2, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                      postHocGroupPerSesNoAdj <- (emmeans (LME_variable2, list (pairwise ~ genotype | session_nr), adjust = "no")) # posthoc test: group per session effect with p value
                      postHocGroupHolmAdj <- emmeans (LME_variable2, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                      postHocGroupPerSesHolmAdj <- (emmeans (LME_variable2, list (pairwise ~ genotype | session_nr), adjust = "holm")) # posthoc test: group per session effect with p value
                      postHocGroupFDRAdj <- emmeans (LME_variable2, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                      postHocGroupPerSesFDRAdj <- (emmeans (LME_variable2, list (pairwise ~ genotype | session_nr), adjust = "FDR")) # posthoc test: group per session effect with p value
                      
                      # combs <- combn(4, 2)
                      # p_values <- numeric(ncol(combs))
                      # for (j in 1:ncol(combs)) {
                      #   if (combs[1, j] == 1) {
                      #     nams1 <- (names(fixef(LME_variable2)))
                      #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                      #     nam_group <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                      #     p_values[j] <- anova(LME_variable2, L = L)$'p-value'      
                      #   } else {
                      #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                      #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                      #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                      #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                      #     p_values[j] <- anova(LME_variable2, L = L)$'p-value'      
                      #   }
                      # }
                      # postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                      # 
                      # export all output of LME_variable2 to textfile
                      outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupPerSesNoAdj, postHocGroupHolmAdj, postHocGroupPerSesHolmAdj, postHocGroupFDRAdj, postHocGroupPerSesFDRAdj)
                      cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model2", ".txt")), sep = "\n", append = FALSE)
                    }
                    
                  
                  # adjusted p-values. multiple comparisons adjustment
                  p.values <- c(0.5048, 0.7087, 0.5439, 0.6799, 0.4459, 0.6911, 0.0528, 0.9702, 0.3153, 0.0660) # outcome measures of unadjusted p.values
                  p.adjust(p.values, method = "holm")
                  
                  
                  
                    # plot model 2 and export
                    # { 
                    # # effect plot
                    #   if (length (unique (ebc$session_nr2)) >= 1) { # effect plot only shown if there is more than 1 session.
                    #     effectPlotData <- function (object, newdata, orig_data) {
                    #       form <- formula(object)
                    #       namesVars <- all.vars(form)
                    #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #       Terms <- delete.response(terms(form))
                    #       mfX <- model.frame(Terms, data = orig_data)
                    #       Terms_new <- attr(mfX, "terms")
                    #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #       X <- model.matrix(Terms_new, mfX_new)
                    #       pred <- c(X %*% betas)
                    #       ses <- sqrt(diag(X %*% V %*% t(X)))
                    #       newdata$pred <- pred
                    #       newdata$low <- pred - 1.96 * ses
                    #       newdata$upp <- pred + 1.96 * ses
                    #       newdata
                    #     }
                    #     newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr)))
                    #     plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                    #                     data = effectPlotData(LME_variable2, newDF, ebc),
                    #                     lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                    #                     xlab = "Session",
                    #                     ylab = (paste0 (variable)))
                    #     print(plot1)
                    #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model2", ".pdf")))
                    #     while (dev.cur()>1) dev.off()
                    #     if(!is.null(dev.list())) dev.off()
                    #   } else {
                    #   }
                    #   # Residual plot
                    #   rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                    #     form <- formula(object)
                    #     namesVars <- all.vars(form)
                    #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #     Terms <- delete.response(terms(form))
                    #     mfX <- model.frame(Terms, data = orig_data)
                    #     Terms_new <- attr(mfX, "terms")
                    #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #     X <- model.matrix(Terms_new, mfX_new)
                    #     ind <- combn(nrow(X), 2)
                    #     k <- ncol(ind)
                    #     out <- matrix(0, k, 5)
                    #     for (i in seq_len(k)) {
                    #       XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                    #       est <- drop(XX %*% betas)
                    #       se <- sqrt(diag(XX %*% V %*% t(XX)))
                    #       out[i, 1] <- est
                    #       out[i, 2] <- se
                    #       out[i, 3] <- est - 1.96 * se
                    #       out[i, 4] <- est + 1.96 * se
                    #       out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                    #     }
                    #     if (k > 2 && adjust.p) {
                    #       out[, 5] <- p.adjust(out[, 5], ...)
                    #     }
                    #     colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                    #     rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                    #     out
                    #   }
                    # 
                    #   v_pVal <- vector(mode="numeric", length=0)
                    #   v_outRowDiff <- vector(mode="numeric", length=0)
                    # 
                    #   for (k in 1:max(ebc$session_nr2))
                    #   {
                    #     # test the effect at specific sessions
                    #     nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype)))
                    #     outRowDiff <- capture.output(nDF, rowDiff(LME_variable2, nDF, ebc))
                    #     v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                    #     pVal <- rowDiff(LME_variable2, nDF, ebc) [5]
                    #     v_pVal <- append (v_pVal, pVal)
                    #   }
                    # 
                    #   v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                    # 
                    #   cat (paste0("Title:", variable, "_session_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_session_diff", "_model2", ".txt")), sep = "\n", append = FALSE)
                    # 
                    #   # plot residuals
                    #   plot2 <- plot(LME_variable2)
                    #   print(plot2)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model2", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # 
                    #   plot3 <- qqnorm(LME_variable2)
                    #   print(plot3)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model2", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # }
                  
                  if ((grepl("time", variable) == TRUE) | (grepl("onset", variable) == TRUE)) {
                  LME_variable2c <-  lme (formula (paste ((colnames(ebc)[i]), "~ genotype + session_nr + session_nr:genotype + amplitude_covariate")),
                                    control = ctrl, data = ebc, random = list(mouse_id = pdDiag(form = ~ poly(session_nr2, 2))), method = "REML", na.action=na.exclude)
                
                    # show results of model 2c and export
                    {
                      sum <- (summary (LME_variable2c))
                      an <- (anova (LME_variable2c))
                      vc <- (VarCorr (LME_variable2c))
                      postHocGroupNoAdj <- emmeans (LME_variable2c, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                      postHocGroupPerSesNoAdj <- (emmeans (LME_variable2c, list (pairwise ~ genotype | session_nr), adjust = "no")) # posthoc test: group per session effect with p value
                      postHocGroupHolmAdj <- emmeans (LME_variable2c, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                      postHocGroupPerSesHolmAdj <- (emmeans (LME_variable2c, list (pairwise ~ genotype | session_nr), adjust = "holm")) # posthoc test: group per session effect with p value
                      postHocGroupFDRAdj <- emmeans (LME_variable2c, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                      postHocGroupPerSesFDRAdj <- (emmeans (LME_variable2c, list (pairwise ~ genotype | session_nr), adjust = "FDR")) # posthoc test: group per session effect with p value
                      
                      # combs <- combn(4, 2)
                      # p_values <- numeric(ncol(combs))
                      # for (j in 1:ncol(combs)) {
                      #   if (combs[1, j] == 1) {
                      #     nams1 <- (names(fixef(LME_variable2c)))
                      #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                      #     nam_group <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                      #     p_values[j] <- anova(LME_variable2c, L = L)$'p-value'      
                      #   } else {
                      #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                      #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                      #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                      #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                      #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                      #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                      #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                      #     p_values[j] <- anova(LME_variable2c, L = L)$'p-value'      
                      #   }
                      # }
                      # postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                      # 
                      # export all output of LME_variable2c to textfile
                      outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupPerSesNoAdj, postHocGroupHolmAdj, postHocGroupPerSesHolmAdj, postHocGroupFDRAdj, postHocGroupPerSesFDRAdj)
                      cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model2c", ".txt")), sep = "\n", append = FALSE)
                    }
                    
                    # plot model 2c and export
                    # { 
                    # # effect plot
                    #   if (length (unique (ebc$session_nr2)) >= 1) { # effect plot only shown if there is more than 1 session.
                    #     effectPlotData <- function (object, newdata, orig_data) {
                    #       form <- formula(object)
                    #       namesVars <- all.vars(form)
                    #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #       Terms <- delete.response(terms(form))
                    #       mfX <- model.frame(Terms, data = orig_data)
                    #       Terms_new <- attr(mfX, "terms")
                    #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #       X <- model.matrix(Terms_new, mfX_new)
                    #       pred <- c(X %*% betas)
                    #       ses <- sqrt(diag(X %*% V %*% t(X)))
                    #       newdata$pred <- pred
                    #       newdata$low <- pred - 1.96 * ses
                    #       newdata$upp <- pred + 1.96 * ses
                    #       newdata
                    #     }
                    #     newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr), amplitude_covariate = median (amplitude_covariate) ))
                    #     plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                    #                     data = effectPlotData(LME_variable2c, newDF, ebc),
                    #                     lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                    #                     xlab = "Session",
                    #                     ylab = (paste0 (variable)))
                    #     print(plot1)
                    #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model2c", ".pdf")))
                    #     while (dev.cur()>1) dev.off()
                    #     if(!is.null(dev.list())) dev.off()
                    #   } else {
                    # 
                    #   }
                    # 
                    #   # Residual plot
                    #   rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                    #     form <- formula(object)
                    #     namesVars <- all.vars(form)
                    #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                    #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                    #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                    #     Terms <- delete.response(terms(form))
                    #     mfX <- model.frame(Terms, data = orig_data)
                    #     Terms_new <- attr(mfX, "terms")
                    #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                    #     X <- model.matrix(Terms_new, mfX_new)
                    #     ind <- combn(nrow(X), 2)
                    #     k <- ncol(ind)
                    #     out <- matrix(0, k, 5)
                    #     for (i in seq_len(k)) {
                    #       XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                    #       est <- drop(XX %*% betas)
                    #       se <- sqrt(diag(XX %*% V %*% t(XX)))
                    #       out[i, 1] <- est
                    #       out[i, 2] <- se
                    #       out[i, 3] <- est - 1.96 * se
                    #       out[i, 4] <- est + 1.96 * se
                    #       out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                    #     }
                    #     if (k > 2 && adjust.p) {
                    #       out[, 5] <- p.adjust(out[, 5], ...)
                    #     }
                    #     colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                    #     rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                    #     out
                    #   }
                    # 
                    #   v_pVal <- vector(mode="numeric", length=0)
                    #   v_outRowDiff <- vector(mode="numeric", length=0)
                    # 
                    #   for (k in 1:max(ebc$session_nr2))
                    #   {
                    #     # test the effect at specific sessions
                    #     nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype), amplitude_covariate = median (amplitude_covariate)))
                    #     outRowDiff <- capture.output(nDF, rowDiff(LME_variable2c, nDF, ebc))
                    #     v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                    #     pVal <- rowDiff(LME_variable2c, nDF, ebc) [5]
                    #     v_pVal <- append (v_pVal, pVal)
                    #   }
                    # 
                    #   v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                    # 
                    #   cat (paste0("Title:", variable, "_session_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_session_diff", "_model2c", ".txt")), sep = "\n", append = FALSE)
                    # 
                    #   # plot residuals
                    #   plot2 <- plot(LME_variable2c)
                    #   print(plot2)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model2c", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # 
                    #   plot3 <- qqnorm(LME_variable2c)
                    #   print(plot3)
                    #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model2c", ".pdf")))
                    #   while (dev.cur()>1) dev.off()
                    #   if(!is.null(dev.list())) dev.off()
                    # } #disabled
                  
                  } else {}
            
                 
                # Check which model fits best and export:
                {      
                  anova(LME_variable0, LME_variable1)
                  anova(LME_variable0, LME_variable2)
                  anova(LME_variable1, LME_variable2)
                  # export all output to textfile
                  ModelCompar <- capture.output(anova(LME_variable0, LME_variable1), anova(LME_variable0, LME_variable2), anova(LME_variable1, LME_variable2))
                  cat (paste0("Title:", variable, "_lme"), ModelCompar, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_ModelCompar", ".txt")), sep = "\n", append = FALSE)
                }

              }
              
              # if there is only 1 session to be analyzed, use simple model 3
              else if (length (unique (ebc$session_nr2)) == 1) 
                {
              
                # model 3: compound symmetry, 1 session.
                LME_variable3 <-  lme (formula (paste ((colnames(ebc)[i]), "~ genotype")),
                                       control = ctrl, data = ebc, correlation = NULL, random = ~ 1 | mouse_id, method = "REML", na.action=na.exclude)
                
                # show results of model 3 and export
                {
                  sum <- (summary (LME_variable3))
                  an <- (anova (LME_variable3))
                  vc <- (VarCorr (LME_variable3))
                  postHocGroupNoAdj <- emmeans (LME_variable3, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                  postHocGroupHolmAdj <- emmeans (LME_variable3, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                  postHocGroupFDRAdj <- emmeans (LME_variable3, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                  
                  # combs <- combn(4, 2)
                  # p_values <- numeric(ncol(combs))
                  # for (j in 1:ncol(combs)) {
                  #   if (combs[1, j] == 1) {
                  #     nams1 <- (names(fixef(LME_variable3)))
                  #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                  #     nam_group <- paste0("genotype", combs[2, j], ":")
                  #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                  #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                  #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                  #     p_values[j] <- anova(LME_variable3, L = L)$'p-value'      
                  #   } else {
                  #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                  #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                  #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                  #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                  #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                  #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                  #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                  #     p_values[j] <- anova(LME_variable3, L = L)$'p-value'      
                  #   }
                  # }
                  # postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                  # 
                  # export all output of LME_variable3 to textfile
                  outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupHolmAdj, postHocGroupFDRAdj)
                  cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model3", ".txt")), sep = "\n", append = FALSE)
                }
                
                # plot model 3 and export
                # { # effect plot
                #   if (length (unique (ebc$session_nr2)) > 1) { # effect plot only shown if there is more than 1 session.
                #     effectPlotData <- function (object, newdata, orig_data) {
                #       form <- formula(object)
                #       namesVars <- all.vars(form)
                #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                #       Terms <- delete.response(terms(form))
                #       mfX <- model.frame(Terms, data = orig_data)
                #       Terms_new <- attr(mfX, "terms")
                #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                #       X <- model.matrix(Terms_new, mfX_new)
                #       pred <- c(X %*% betas)
                #       ses <- sqrt(diag(X %*% V %*% t(X)))
                #       newdata$pred <- pred
                #       newdata$low <- pred - 1.96 * ses
                #       newdata$upp <- pred + 1.96 * ses
                #       newdata
                #     }
                #     newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr)))
                #     plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                #                     data = effectPlotData(LME_variable3, newDF, ebc),
                #                     lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                #                     xlab = "Session",
                #                     ylab = (paste0 (variable)))
                #     print(plot1)
                #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model3", ".pdf")))
                #     while (dev.cur()>1) dev.off()
                #     if(!is.null(dev.list())) dev.off()
                #   } else {
                #   }
                #   # Residual plot 
                #   rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                #     form <- formula(object)
                #     namesVars <- all.vars(form)
                #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                #     Terms <- delete.response(terms(form))
                #     mfX <- model.frame(Terms, data = orig_data)
                #     Terms_new <- attr(mfX, "terms")
                #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                #     X <- model.matrix(Terms_new, mfX_new)
                #     ind <- combn(nrow(X), 2)
                #     k <- ncol(ind)
                #     out <- matrix(0, k, 5)
                #     for (i in seq_len(k)) {
                #       XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                #       est <- drop(XX %*% betas)
                #       se <- sqrt(diag(XX %*% V %*% t(XX)))
                #       out[i, 1] <- est
                #       out[i, 2] <- se
                #       out[i, 3] <- est - 1.96 * se
                #       out[i, 4] <- est + 1.96 * se
                #       out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                #     }
                #     if (k > 2 && adjust.p) {
                #       out[, 5] <- p.adjust(out[, 5], ...)
                #     }
                #     colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                #     rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                #     out
                #   }
                #   
                #   v_pVal <- vector(mode="numeric", length=0)
                #   v_outRowDiff <- vector(mode="numeric", length=0)
                #   
                #   for (k in 1:max(ebc$session_nr2)) 
                #   {
                #     # test the effect at specific sessions
                #     nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype)))
                #     outRowDiff <- capture.output(nDF, rowDiff(LME_variable3, nDF, ebc))
                #     v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                #     pVal <- rowDiff(LME_variable3, nDF, ebc) [5]
                #     v_pVal <- append (v_pVal, pVal)
                #   }
                #   
                #   v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                #   
                #   cat (paste0("Title:", variable, "_group_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_group_diff", "_model3", ".txt")), sep = "\n", append = FALSE)
                #   
                #   # plot residuals
                #   plot2 <- plot(LME_variable3)
                #   print(plot2)
                #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model3", ".pdf")))
                #   while (dev.cur()>1) dev.off()
                #   if(!is.null(dev.list())) dev.off()
                # 
                #   plot3 <- qqnorm(LME_variable3)
                #   print(plot3)
                #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model3", ".pdf")))
                #   while (dev.cur()>1) dev.off()
                #   if(!is.null(dev.list())) dev.off()
                # } 
              
              
              # if 'time' is in columnheader, test same model 0 but now use AmpattX as covariate.  
              if ((grepl("time", variable) == TRUE) | (grepl("onset", variable) == TRUE)) {
                # Model 0c: idem as model 0 but now with AmpattX as fixed factor added
                LME_variable3c <- lme (formula (paste ((colnames(ebc)[i]), "~ genotype + amplitude_covariate")),
                                       control = ctrl, data = ebc, correlation = NULL, random = ~ 1 | mouse_id, method = "REML", na.action=na.exclude) 
                
                # show results of model 3c and export
                {
                  sum <- (summary (LME_variable3c))
                  an <- (anova (LME_variable3c))
                  vc <- (VarCorr (LME_variable3c))
                  postHocGroupNoAdj <- emmeans (LME_variable3c, list (pairwise ~ genotype), adjust = "no") # posthoc test: overall group effect with p value
                  postHocGroupHolmAdj <- emmeans (LME_variable3c, list (pairwise ~ genotype), adjust = "holm") # posthoc test: overall group effect with p value
                  postHocGroupFDRAdj <- emmeans (LME_variable3, list (pairwise ~ genotype), adjust = "FDR") # posthoc test: overall group effect with p value
                 
                  #  combs <- combn(4, 2)
                  # p_values <- numeric(ncol(combs))
                  # for (j in 1:ncol(combs)) {
                  #   if (combs[1, j] == 1) {
                  #     nams1 <- (names(fixef(LME_variable3c)))
                  #     nams <- (nams1[!nams1 %in% "amplitude_covariate"])
                  #     nam_group <- paste0("genotype", combs[2, j], ":")
                  #     grp_ind <- grep(nam_group, nams, fixed = TRUE)
                  #     L <- matrix(0.0, nrow = length(grp_ind), ncol = length(nams))
                  #     L[cbind(1:length(grp_ind), grp_ind)] <- 1 
                  #     p_values[j] <- anova(LME_variable3c, L = L)$'p-value'      
                  #   } else {
                  #     nam_group_i <- paste0("genotype", combs[1, j], ":")
                  #     grp_ind_i <- grep(nam_group_i, nams, fixed = TRUE)
                  #     nam_group_j <- paste0("genotype", combs[2, j], ":")
                  #     grp_ind_j <- grep(nam_group_j, nams, fixed = TRUE)
                  #     L <- matrix(0.0, nrow = length(grp_ind_i), ncol = length(nams))
                  #     L[cbind(1:length(grp_ind_i), grp_ind_i)] <- 1
                  #     L[cbind(1:length(grp_ind_j), grp_ind_j)] <- -1
                  #     p_values[j] <- anova(LME_variable3c, L = L)$'p-value'      
                  #   }
                  # }
                  # postHocgenotypees <- round(cbind(t(combs), p_values, p.adjust(p_values, method = "holm"), p.adjust(p_values, method = "fdr")),  4)
                  # 
                  # export all output of LME_variable3c to textfile
                  outLME <- capture.output(sum, an, vc, postHocGroupNoAdj, postHocGroupHolmAdj, postHocGroupFDRAdj)
                  cat (paste0("Title:", variable, "_lme"), outLME, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_model3c", ".txt")), sep = "\n", append = FALSE)
                }
                
                # plot model 3c and export
                # { # effect plot
                #   if (length (unique (ebc$session_nr2)) > 1) { # effect plot only shown if there is more than 1 session.
                #     effectPlotData <- function (object, newdata, orig_data) {
                #       form <- formula(object)
                #       namesVars <- all.vars(form)
                #       betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                #       V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                #       orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                #       Terms <- delete.response(terms(form))
                #       mfX <- model.frame(Terms, data = orig_data)
                #       Terms_new <- attr(mfX, "terms")
                #       mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                #       X <- model.matrix(Terms_new, mfX_new)
                #       pred <- c(X %*% betas)
                #       ses <- sqrt(diag(X %*% V %*% t(X)))
                #       newdata$pred <- pred
                #       newdata$low <- pred - 1.96 * ses
                #       newdata$upp <- pred + 1.96 * ses
                #       newdata
                #     }
                #     newDF <- with(ebc, expand.grid(genotype = levels(genotype),  session_nr = levels(session_nr), amplitude_covariate = median (amplitude_covariate) ))
                #     plot1 <- xyplot(pred + low + upp ~ session_nr | genotype,
                #                     data = effectPlotData(LME_variable3c, newDF, ebc),
                #                     lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                #                     xlab = "Session",
                #                     ylab = (paste0 (variable)))
                #     print(plot1)
                #     dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_xyplot_95CI_", "model3c", ".pdf")))
                #     while (dev.cur()>1) dev.off()
                #     if(!is.null(dev.list())) dev.off()
                #   } else {
                #     
                #   }
                #   
                #   # Residual plot 
                #   rowDiff <- function (object, newdata, orig_data, adjust.p = FALSE, ...) {
                #     form <- formula(object)
                #     namesVars <- all.vars(form)
                #     betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
                #     V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
                #     orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
                #     Terms <- delete.response(terms(form))
                #     mfX <- model.frame(Terms, data = orig_data)
                #     Terms_new <- attr(mfX, "terms")
                #     mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
                #     X <- model.matrix(Terms_new, mfX_new)
                #     ind <- combn(nrow(X), 2)
                #     k <- ncol(ind)
                #     out <- matrix(0, k, 5)
                #     for (i in seq_len(k)) {
                #       XX <- X[ind[1, i], , drop = FALSE] - X[ind[2, i], , drop = FALSE]
                #       est <- drop(XX %*% betas)
                #       se <- sqrt(diag(XX %*% V %*% t(XX)))
                #       out[i, 1] <- est
                #       out[i, 2] <- se
                #       out[i, 3] <- est - 1.96 * se
                #       out[i, 4] <- est + 1.96 * se
                #       out[i, 5] <- 2 * pnorm(abs(est / se), lower.tail = FALSE)
                #     }
                #     if (k > 2 && adjust.p) {
                #       out[, 5] <- p.adjust(out[, 5], ...)
                #     }
                #     colnames(out) <- c("Diff", "Std.Err.", "95%low", "95%upp", "p-value")
                #     rownames(out) <- paste(ind[1, ], "-", ind[2, ])
                #     out
                #   }
                #   
                #   v_pVal <- vector(mode="numeric", length=0)
                #   v_outRowDiff <- vector(mode="numeric", length=0)
                #   
                #   for (k in 1:max(ebc$session_nr2)) 
                #   {
                #     # test the effect at specific sessions
                #     nDF <- with(ebc, expand.grid(session_nr = factor(k, levels(session_nr)), genotype = levels(genotype), amplitude_covariate = median (amplitude_covariate) ))
                #     outRowDiff <- capture.output(nDF, rowDiff(LME_variable3c, nDF, ebc))
                #     v_outRowDiff <- append (v_outRowDiff, outRowDiff)
                #     pVal <- rowDiff(LME_variable3c, nDF, ebc) [5]
                #     v_pVal <- append (v_pVal, pVal)
                #   }
                #   
                #   v_pValAdj <- p.adjust(v_pVal) # adjust for multiple comparison using default Holms correction
                #   
                #   cat (paste0("Title:", variable, "_group_diff"), "\n", v_outRowDiff, "\n", "Holms corrected p-values", v_pValAdj, file = (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_group_diff", "_model3c", ".txt")), sep = "\n", append = FALSE)
                #   
                #   # plot residuals
                #   plot2 <- plot(LME_variable3c)
                #   print(plot2)
                #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_residualPlot", "_model3c", ".pdf")))
                #   while (dev.cur()>1) dev.off()
                #   if(!is.null(dev.list())) dev.off()
                # 
                #   plot3 <- qqnorm(LME_variable3c)
                #   print(plot3)
                #   dev.copy (pdf, (paste0 (wd, "/", foldername, "/", variable, "_", i, "/", variable, "_qqnorm_residuals", "_model3c", ".pdf")))
                #   while (dev.cur()>1) dev.off()
                #   if(!is.null(dev.list())) dev.off()
                # } 
              
                } else {
                
              }

              } 
              
              else {} #,main = i, error=function(e) {print ("unspecified error")})
      }
    }, main = i, error=function(e) {print ("unspecified error")})
  }
}


################
##AD HOC PLOTS##
################

# make new folder in the wd and use it to store output
dir.create(paste0 (wd, "/R-analysis graphics"), showWarnings = FALSE)

# scatter1 ISIpeakamp - ISIpeaktime ~ group
{scatter1 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
  geom_point (aes(color=genotype)) +
  geom_smooth(method=loess, se = TRUE) +
  scale_x_continuous(limits = c(500,1000)) +
  scale_y_continuous(limits = c(0,2)) +
  facet_grid(.~ genotype) +
  theme_classic()

print(scatter1)
dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatter1", ".pdf")))
while (dev.cur()>1) dev.off()
if(!is.null(dev.list())) dev.off()}

# scatter2 ISIpeakamp - ISIpeaktime ~ group + session
{scatter2 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    geom_point (aes(color=genotype)) +
    geom_smooth(method=loess, se = TRUE) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(scatter2)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatter2", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatter3 CRpeakamp - CRpeaktime ~ group
{scatter3 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_point (aes(color=genotype)) +
    geom_smooth(method=loess, se = TRUE) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(.~ genotype) +
    theme_classic()
  
  print(scatter3)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatter3", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatter4 CRpeakamp - CRpeaktime ~ group + session
{scatter4 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_point (aes(color=genotype)) +
    geom_smooth(method=loess, se = TRUE) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(scatter4)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatter4", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterheat1 ISIpeakamp - ISIpeaktime ~ group
{scatterheat1 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    geom_point (color = "black") +
    stat_density_2d(aes(fill = ..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="red") +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid(.~genotype) +
    theme_classic()
  
  print(scatterheat1)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterheat1", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterheat2 ISIpeakamp - ISIpeaktime ~ group + session
{scatterheat2 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    geom_point (color = "black") +
    stat_density_2d(aes(fill = ..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="red") +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid (session_nr ~ genotype) +
    theme_classic()
  
  print(scatterheat2)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterheat2", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterheat3 CRpeakamp - CRpeaktime ~ group
{scatterheat3 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_point (color = "black") +
    stat_density_2d(aes(fill = ..level..), geom="polygon") +
    scale_fill_gradient(low="blue", high="red") +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid(.~genotype) +
    theme_classic()
  
  print(scatterheat3)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterheat3", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterheat4 CRpeakamp - CRpeaktime ~ group + session
{scatterheat4 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_point (color = "black") +
    stat_density_2d (geom="polygon") +
    scale_fill_gradientn (colours=jet.colors(7), legend_param=list(colorbar=T, colorbar_nbin=100)) +
    stat_smooth() +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid (.~ genotype) +
    theme_classic()
  
  print(scatterheat4)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterheat4", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterline1 ISIpeakamp - ISIpeaktime ~ group 
{scatterline1 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    geom_point (color = "black") +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(.~ genotype) +
    theme_classic()
  
  print(scatterline1)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterline1", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterline2 ISIpeakamp - ISIpeaktime ~ group + session
{scatterline2 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    geom_point (color = "black") +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(scatterline2)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterline2", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterline3 CRpeakamp - CRpeaktime ~ group 
{scatterline3 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_point (color = "black") +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(.~ genotype) +
    theme_classic()
  
  print(scatterline3)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterline3", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# scatterline4 CRpeakamp - CRpeaktime ~ group + session
{
  scatterline4 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_point (color = "black") +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(scatterline4)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "scatterline4", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}


# lineline1 ISIpeakamp - ISIpeaktime ~ group 
{lineline1 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(.~ genotype) +
    theme_classic()
  
  print(lineline1)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "lineline1", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# lineline2 ISIpeakamp - ISIpeaktime ~ group + session
{lineline2 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(lineline2)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "lineline2", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# lineline3 CRpeakamp - CRpeaktime ~ group 
{lineline3 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(.~ genotype) +
    theme_classic()
  
  print(lineline3)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "lineline3", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# lineline4 CRpeakamp - CRpeaktime ~ group + session
{
  lineline4 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    geom_density_2d (alpha = 1) +
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(0,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(lineline4)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "lineline4", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}


jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))

# heat1 ISIpeakamp - ISIpeaktime ~ group
{heat1 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
    scale_fill_gradientn(colours=jet.colors(9)) + 
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid(.~genotype) +
    theme_classic()
  
  print(heat1)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "heat1", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# heat2 ISIpeakamp - ISIpeaktime ~ group
{heat2 <- ggplot(ebcAll, aes(x=ISIpeaktime, y=ISIpeakamp, group=genotype)) +
    stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
    scale_fill_gradientn(colours=jet.colors(9)) + 
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(heat2)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "heat2", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# heat3 CRpeakamp - CRpeaktime ~ group
{heat3 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
    scale_fill_gradientn(colours=jet.colors(9)) + 
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid(.~genotype) +
    theme_classic()
  
  print(heat3)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "heat3", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}

# heat4 CRpeakamp - CRpeaktime ~ group
{heat4 <- ggplot(ebcAll, aes(x=CRpeaktime, y=CRpeakamp, group=genotype)) +
    stat_density2d(geom="tile", aes(fill = ..density..), contour = FALSE) +
    scale_fill_gradientn(colours=jet.colors(9)) + 
    scale_x_continuous(limits = c(500,1000)) +
    scale_y_continuous(limits = c(-.25,2)) +
    facet_grid(session_nr ~ genotype) +
    theme_classic()
  
  print(heat4)
  dev.copy (pdf, (paste0 (wd, "/R-analysis graphics/", "heat4", ".pdf")))
  while (dev.cur()>1) dev.off()
  if(!is.null(dev.list())) dev.off()}
 









