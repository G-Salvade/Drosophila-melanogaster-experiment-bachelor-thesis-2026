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
WMOA <- read.csv2(file.choose())

# ANOVA and assumption
M_WMOA <- lm(Average ~ Karyotype, WMOA)
anova(M_WMOA)
shapiro.test(aov(M_WMOA)$residuals)
WMOA$Average <- as.numeric(as.character(WMOA$Average))
str(WMOA)
leveneTest(Average ~ Karyotype, WMOA)

# Histogramm and QQ Plot
model <- aov(M_WMOA)
par(mar = c(5, 5, 4, 2))
qqPlot(aov(M_WMOA)$residuals,
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


#preparing differents Karyotype
WMOA_3RP <- WMOA[c(which(str_detect(WMOA$Karyotype, "3RP|STD"))),]
WMOA_2Lt <- WMOA[c(which(str_detect(WMOA$Karyotype, "2Lt|STD"))),]
WMOA_3RK <- WMOA[c(which(str_detect(WMOA$Karyotype, "3RK|STD"))),]
WMOA_HOMO <- WMOA[-c(which(str_detect(WMOA$Karyotype, "HET"))),]


# Box Plot
WMOA_B1 <- ggplot(WMOA_3RP, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.65,0.85) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male wing area 3RP-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

WMOA_B2 <- ggplot(WMOA_3RK, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.65,0.85) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male wing area 3RK-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


WMOA_B3 <- ggplot(WMOA_2Lt, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.65,0.85) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male wing area 2Lt-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


WMOA_B4 <- ggplot(WMOA_HOMO, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.65,0.85) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male wing area (Homokaryotype)",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Wing area (") * bold(mm^2) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

# General box plot
(WMOA_B1|WMOA_B2)/(WMOA_B3|WMOA_B4) + plot_annotation(title="Male Wing Area", 
                                                      tag_levels='A', 
                                                      theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12))) 
                                                      



