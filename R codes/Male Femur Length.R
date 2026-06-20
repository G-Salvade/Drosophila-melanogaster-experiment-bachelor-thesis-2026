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
MFOA <- read.csv2(file.choose())

# ANOVA and assumption
M_MFOA <- lm(Average ~ Karyotype, MFOA)
anova(M_MFOA)

MFOA$Average <- as.numeric(as.character(MFOA$Average))
str(MFOA)
shapiro.test(aov(M_MFOA)$residuals)
leveneTest(Average ~ Karyotype, MFOA)


# Removing outlier
WFOA1 <- WFOA[(WFOA$Average <1.11),]

# ANOVA and assumption
M_WFOA1 <- lm(Average ~ Karyotype, WFOA1)
anova(M_WFOA1)
shapiro.test(aov(M_WFOA1)$residuals)
leveneTest(Average~ Karyotype, WFOA1)

# Histogramm and QQ Plot
model <- aov(M_MFOA)
par(mar = c(5, 5, 4, 2))
qqPlot(aov(M_MFOA)$residuals,
  main = "Normal Q-Q Plot",
  xlab = "Theoretical Quantiles",
  ylab = expression("Residuals Value"),
  id= FALSE) + mtext("A", side = 3, line = 1.4, adj = -0.05,
          font = 2, cex = 1.4)

hist(residuals(model),
  main = "Histogramm",
  xlab = expression("Residuals Value"),
  ylab = "Frequency")mtext("B", side = 3, line = 1.4, adj = -0.05,
      font = 2, cex = 1.4)



# Group Dataset
MFOA_3RP <- MFOA[c(which(str_detect(MFOA$Karyotype, "3RP|STD"))),]
MFOA_2Lt <- MFOA[c(which(str_detect(MFOA$Karyotype, "2Lt|STD"))),]
MFOA_3RK <- MFOA[c(which(str_detect(MFOA$Karyotype, "3RK|STD"))),]
MFOA_HOMO <- MFOA[-c(which(str_detect(MFOA$Karyotype, "HET"))),]


# Box Plot
MFOA_B1 <- ggplot(MFOA_3RP, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.5,0.58) +stat_pvalue_manual(p_MFOA_3RP, label="p.signif", hide.ns = TRUE, y.position=0.57, step.increase = 0.01)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male femur length 3RP-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

MFOA_B2 <- ggplot(MFOA_3RK, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.5,0.58) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male femur length 3RK-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

MFOA_B3 <- ggplot(MFOA_2Lt, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.5,0.58)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male femur length 2Lt-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


MFOA_B4 <- ggplot(MFOA_HOMO, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.5,0.58) +stat_pvalue_manual(p_MFOA_HOMO, label="p.signif", hide.ns = TRUE, y.position=0.57, step.increase = 0.01)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male femur length (Homokaryotype)",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

# General box plot
(MFOA_B1|MFOA_B2)/(MFOA_B3|MFOA_B4) + plot_annotation(title="Male Femur Length", 
                                                      tag_levels='A', 
                                                      theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))


                                                  

