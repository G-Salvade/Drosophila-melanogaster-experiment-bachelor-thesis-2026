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
MTOA <- read.csv2(file.choose())

# ANOVA and assumption
M_MTOA <- lm(Average ~ Karyotype, MTOA)
anova(M_MTOA)

MTOA$Average <- as.numeric(as.character(MTOA$Average))
str(MTOA)
shapiro.test(aov(M_MTOA)$residuals)
leveneTest(Average ~ Karyotype, MTOA)

# Histogramm and QQ Plot
model <- aov(M_MTOA)
par(mfrow = c(1, 2))
par(mar = c(5, 5, 4, 2))
qqPlot(aov(M_MTOA)$residuals,
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


# Contrast test: nothing significant
a.18 = 1-(1-0.05)^(1/18)

All.comparisons.3RP.MTOA <- dunn_test(MTOA_3RP, Average ~ Karyotype, p.adjust.method = "none")
All.comparisons.2Lt.MTOA <- dunn_test(MTOA_2Lt, Average ~ Karyotype, p.adjust.method = "none")
All.comparisons.3RK.MTOA <- dunn_test(MTOA_3RK, Average ~ Karyotype, p.adjust.method = "none")

Contrast.MTOA <- rbind(All.comparisons.3RP.MTOA, All.comparisons.2Lt.MTOA, All.comparisons.3RK.MTOA)
S_Contrast.MTOA <- Contrast.MTOA[Contrast.MTOA$p< a.18,]
S_Contrast.MTOA <- mutate(Contrast.MTOA, p.signif =case_when(p< a.18 ~ "*", TRUE ~ "ns"))

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
MTOA_B1 <- ggplot(MTOA_3RP, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.44,0.51) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male tibia length 3RP-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

MTOA_B2 <- ggplot(MTOA_3RK, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.44,0.51) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male tibia length 3RK-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


MTOA_B3 <- ggplot(MTOA_2Lt, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.44,0.51) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male tibia length 2Lt-STD",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


MTOA_B4 <- ggplot(MTOA_HOMO, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.44,0.51) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Male tibia length (Homokaryotype)",subtitle = "\u2642", x = "Karyotype", y = expression(bold("Tibia length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))


# General Box Plot
(MTOA_B1|MTOA_B2)/(MTOA_B3|MTOA_B4) + plot_annotation(title="Male Tibia Length", 
                                                      tag_levels='A', 
                                                      theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))
