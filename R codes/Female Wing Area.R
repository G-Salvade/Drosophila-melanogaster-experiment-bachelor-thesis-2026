library(ggplot2)
library(ggsignif)
library(readxl)
library(FSA)
library(stringr)
library(patchwork)
library(tidyverse)
library(car)
library(lattice)
library(rstatix)
library(ggpubr)
library(stats)
library(emmeans)
library(dplyr)


# Set colors
Colors = c(
  "STD" ="#00b4d8",
  "3RP"="#f25c54",
  "3RP_HET_1"="#ffbf69",
  "3RP_HET_2"="#ffdab9",
  "2Lt"="#52b788",
  "2Lt_HET_1"="#80ed99",
  "2Lt_HET_2"="#c7f9cc",
  "3RK"="#8367c7",
  "3RK_HET_1"="#e0aaff",
  "3RK_HET_2"="#ebd9fc"
)

# Import Data
WFOA <- read.csv2(file.choose())

# ANOVA and assumption
M_WFOA <- lm(Average ~ Karyotype, WFOA)
anova(M_WFOA)

WFOA$Average <- as.numeric(as.character(WFOA$Average))
str(WFOA)
leveneTest(Average ~ Karyotype, WFOA)
shapiro.test(aov(M_WFOA)$residuals)

# Removing outlier
WFOA1 <- WFOA[(WFOA$Average <1.11),]

# ANOVA and assumption
M_WFOA1 <- lm(Average ~ Karyotype, WFOA1)
anova(M_WFOA1)
shapiro.test(aov(M_WFOA1)$residuals)
leveneTest(Average~ Karyotype, WFOA1)

# Histogramm and QQ Plot
model <- aov(M_WFOA1)
par(mar = c(5, 5, 4, 2))
qqPlot(aov(M_WFOA1)$residuals,
  main = "Normal Q-Q Plot",
  xlab = "Theoretical Quantiles",
  ylab = expression("Residuals Value"),
  id= FALSE) + mtext("A", side = 3, line = 1.4, adj = -0.05,
          font = 2, cex = 1.4)

hist(residuals(model),
  main = "Histogramm",
  xlab = expression("Residuals Value"),
  ylab = "Frequency")
mtext("B", side = 3, line = 1.4, adj = -0.05,
      font = 2, cex = 1.4)


# Contrast test for WFOA1
WFOA1_3RP <- WFOA1[c(which(str_detect(WFOA1$Karyotype, "3RP|STD"))),]
WFOA1_2Lt <- WFOA1[c(which(str_detect(WFOA1$Karyotype, "2Lt|STD"))),]
WFOA1_3RK <- WFOA1[c(which(str_detect(WFOA1$Karyotype, "3RK|STD"))),]

M_WFOA1_3RP <- lm(Average ~ Karyotype, WFOA1_3RP)
averageWFOA1_3RP <- emmeans(M_WFOA1_3RP, "Karyotype")
ComparaisonWFOA1_3RP <- pairs(averageWFOA1_3RP, adjust="none")
ComparaisonWFOA1_3RP <-as.data.frame(ComparaisonWFOA1_3RP)

M_WFOA1_2Lt <- lm(Average ~ Karyotype, WFOA1_2Lt)
averageWFOA1_2Lt <- emmeans(M_WFOA1_2Lt, "Karyotype")
ComparaisonWFOA1_2Lt <- pairs(averageWFOA1_2Lt, adjust="none")
ComparaisonWFOA1_2Lt <-as.data.frame(ComparaisonWFOA1_2Lt)

M_WFOA1_3RK <- lm(Average ~ Karyotype, WFOA1_3RK)
averageWFOA1_3RK <- emmeans(M_WFOA1_3RK, "Karyotype")
ComparaisonWFOA1_3RK <- pairs(averageWFOA1_3RK, adjust="none")
ComparaisonWFOA1_3RK <-as.data.frame(ComparaisonWFOA1_3RK)


Total_contrast <- rbind(ComparaisonWFOA1_3RP,ComparaisonWFOA1_2Lt, ComparaisonWFOA1_3RK)
adjusted_p <- p.adjust(Total_contrast$p.value, method = "bonferroni")

# alpha orrection
fwa=1-(1-0.05)^(1/18)

#FINAL RESULTS FOR CONTRAST TEST
Significant_Total_contrast <- Total_contrast[Total_contrast$p.value<fwa,] 

# Add p-value
Contrast_WFOA1<- separate(Total_contrast, contrast, into= c("group1", "group2"), sep = " - ") 
Contrast_WFOA1_sig <- mutate(Contrast_WFOA1, p.signif=case_when(p.value < fwa ~ "*",TRUE ~ "ns"))   

WFOA1_3RP <- WFOA1[c(which(str_detect(WFOA1$Karyotype, "3RP|STD"))),]
WFOA1_2Lt <- WFOA1[c(which(str_detect(WFOA1$Karyotype, "2Lt|STD"))),]
WFOA1_3RK <- WFOA1[c(which(str_detect(WFOA1$Karyotype, "3RK|STD"))),]
WFOA1_HOMO <- WFOA1[-c(which(str_detect(WFOA1$Karyotype, "HET"))),]

p_WFOA1_3RP <- Contrast_WFOA1_sig[Contrast_WFOA1_sig$group1 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2") & Contrast_WFOA1_sig$group2 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2"),]
p_WFOA1_2Lt <- Contrast_WFOA1_sig[Contrast_WFOA1_sig$group1 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2") & Contrast_WFOA1_sig$group2 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2"),]
p_WFOA1_3RK <- Contrast_WFOA1_sig[Contrast_WFOA1_sig$group1 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2") & Contrast_WFOA1_sig$group2 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2"),]
p_WFOA1_HOMO <- Contrast_WFOA1_sig[Contrast_WFOA1_sig$group1 %in% c("STD","3RP","2Lt","3RK") & Contrast_WFOA1_sig$group2 %in% c("STD","3RP","2Lt","3RK"),]


# Box Plot
WFOA1_B1 <- ggplot(WFOA1_3RP, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.8,1.2) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female wing area 3RP-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

WFOA1_B2 <- ggplot(WFOA1_3RK, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.8,1.2) +stat_pvalue_manual(p_WFOA1_3RK, label="p.signif", hide.ns = TRUE, y.position = 1.15, step.increase = 0.08)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female wing area 3RK-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


WFOA1_B3 <- ggplot(WFOA1_2Lt, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.8 ,1.2) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female wing area 2Lt-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

WFOA1_B4 <- ggplot(WFOA1_HOMO, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.8 ,1.2) +stat_pvalue_manual(p_WFOA1_HOMO, label="p.signif", hide.ns = TRUE, y.position = 1.15, step.increase = 0.08)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female wing area (Homokaryotype)",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


# General box plot
(WFOA1_B1|WFOA1_B2)/(WFOA1_B3|WFOA1_B4) + plot_annotation(title="Female Wing Area", 
                                                          tag_levels='A', 
                                                          theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))






