# D409V/WT medial septum MEA analysis 
# Author: Millie Sander
# Part 2: Inter-spike interval 

# loading required packages
library(xlsx)
library(lme4)
library(emmeans)
library(car)
library(fitdistrplus)

ISI <- read.xlsx("/lustre/home/ms739/GBA_MEA/ISI.xlsx", header = TRUE, sheetIndex = "Sheet1")

# creating new dataframe to alter (so dont lose the original)
# changing all NaN values to 1000 (bin size)
data <- ISI
data$isi[data$isi == "NaN"] <- 1000

write.xlsx(data, "/lustre/home/ms739/GBA_MEA/ISI_NaNtransformed.xlsx")

data$Genotype <- as.factor(data$Genotype)
data$Recording <- as.factor(data$Recording)
data$Cell_ID <- as.factor(data$Cell_ID)
data$Animal_ID <- as.factor(data$Animal_ID)
data$Condition <- as.factor(data$Condition)
data$Sex <- as.factor(data$Sex)
data$isi <- as.numeric(data$isi)

# creating linear models with and without sex included as a variable 
m1 <- lmer(isi ~ Genotype*Condition + (1|Animal_ID/Recording/Cell_ID), data = data)
m2 <- lmer(isi ~ Genotype*Condition + Sex + (1|Animal_ID/Recording/Cell_ID), data = data)

# seeing if there is a sig diff between models, and which is the best fit 
anova(m1, m2)

# testing assumptions of chosen model (m1)
qqnorm(residuals(m2))
qqline(residuals(m2))
hist(residuals(m2))
plot(m2)

vif(m2)

leveneTest(residuals(m2)~data$Genotype)
shapiro.test(residuals(m2))

# checking theoretical distribution 
descdist(data$isi)

# beta glmer using scaled Hz variable and glmmTMB package 
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
data$scaled <- scale_values(data$isi)
data$scaled <- as.numeric(data$scaled)
data$scaled[data$scaled == "0"] <- 0.00000001
data$scaled[data$scaled == "1"] <- 0.99

m2b <- glmmTMB(scaled ~ Genotype*Condition + Sex + (1|Animal_ID/Recording/Cell_ID), data = data, family = beta_family())

joint_tests(m2b)

emmeans(m2b, list(pairwise ~ Genotype*Condition), adjust = "bonferroni")

