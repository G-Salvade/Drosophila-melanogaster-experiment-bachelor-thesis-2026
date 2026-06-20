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
FTOA <- read.csv2(file.choose())

# ANOVA and assumption
M_FTOA <- lm(Average ~ Karyotype, FTOA)
anova(M_FTOA)

FTOA$Average <- as.numeric(as.character(FTOA$Average))
str(FTOA)
shapiro.test(aov(M_FTOA)$residuals)
leveneTest(Average ~ Karyotype, FTOA)

# Histogramm and QQ Plot
model <- aov(M_FTOA)
par(mar = c(5, 5, 4, 2))
qqPlot(aov(M_FTOA)$residuals,
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


# Non parametric Test 
kruskal.test(Average ~ Karyotype, FTOA)

# Contrast test 
All.comparisons.3RP.FTOA <- dunn_test(FTOA_3RP, Average ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.FTOA <- dunn_test(FTOA_2Lt, Average ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.FTOA <- dunn_test(FTOA_3RK, Average ~ Karyotype, p.adjust.method = "none")

a.18 = 1-(1-0.05)^(1/18)

Contrast.FTOA <- rbind(All.comparisons.3RP.FTOA, All.comparisons.2Lt.FTOA, All.comparisons.3RK.FTOA)
S_Contrast.FTOA <- Contrast.FTOA[Contrast.FTOA$p< a.18,]
S_Contrast.FTOA <- mutate(Contrast.FTOA, p.signif =case_when(p< a.18 ~ "*", TRUE ~ "ns"))

# Add p-value
FTOA_3RP <- FTOA[c(which(str_detect(FTOA$Karyotype, "3RP|STD"))),]
FTOA_2Lt <- FTOA[c(which(str_detect(FTOA$Karyotype, "2Lt|STD"))),]
FTOA_3RK <- FTOA[c(which(str_detect(FTOA$Karyotype, "3RK|STD"))),]
FTOA_HOMO <- FTOA[-c(which(str_detect(FTOA$Karyotype, "HET"))),]

p_FTOA_3RP <- S_Contrast.FTOA[S_Contrast.FTOA$group1 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2") & S_Contrast.FTOA$group2 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2"),]
p_FTOA_2Lt <- S_Contrast.FTOA[S_Contrast.FTOA$group1 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2") & S_Contrast.FTOA$group2 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2"),]
p_FTOA_3RK <-S_Contrast.FTOA[S_Contrast.FTOA$group1 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2") & S_Contrast.FTOA$group2 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2"),]
p_FTOA_HOMO <- S_Contrast.FTOA[S_Contrast.FTOA$group1 %in% c("STD","3RP","2Lt","3RK") & S_Contrast.FTOA$group2 %in% c("STD","3RP","2Lt","3RK"),]

p_FTOA_3RP_sig <- subset(p_FTOA_3RP, p.signif != "ns")
p_FTOA_3RK_sig <- subset(p_FTOA_3RK, p.signif != "ns")
p_FTOA_2Lt_sig <- subset(p_FTOA_2Lt, p.signif != "ns")
p_FTOA_HOMO_sig <- subset(p_FTOA_HOMO, p.signif != "ns")


# Box Plot
FTOA_B1 <- ggplot(FTOA_3RP, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.455,0.521) +stat_pvalue_manual(p_FTOA_3RP_sig, label="p.signif", hide.ns = TRUE, y.position=0.505, step.increase = 0.08)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female tibia length 3RP-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


FTOA_B2 <- ggplot(FTOA_3RK, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.455,0.521) +stat_pvalue_manual(p_FTOA_3RK_sig, label="p.signif", hide.ns = TRUE, y.position=0.505, step.increase = 0.08)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female tibia length 3RK-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


FTOA_B3 <- ggplot(FTOA_2Lt, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.455,0.521) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female tibia length 2Lt-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


FTOA_B4 <- ggplot(FTOA_HOMO, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.455,0.521) +stat_pvalue_manual(p_FTOA_HOMO_sig, label="p.signif", hide.ns = TRUE, y.position=0.505, step.increase = 0.08)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female tibia length (Homokaryotype)",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))



# General Box Plot
(FTOA_B1|FTOA_B2)/(FTOA_B3|FTOA_B4) + plot_annotation(title="Female Tibia Length", 
                                                      tag_levels='A', 
                                                      theme=theme(plot.title=element_text(hjust=0.48, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))


