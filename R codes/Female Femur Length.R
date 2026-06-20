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
FFOA <- read.csv2(file.choose())

# ANOVA and assumption
M_FFOA <- lm(Average ~ Karyotype, FFOA)
anova(M_FFOA)

FFOA$Average <- as.numeric(as.character(FFOA$Average))
str(FFOA)
shapiro.test(aov(M_FFOA)$residuals)
leveneTest(Average ~ Karyotype, FFOA)

# Removing outlier
FFOA1 <- FFOA[(FFOA$Average <0.608),]

# ANOVA and assumption
M_FFOA1 <- lm(Average ~ Karyotype, FFOA1)
anova(M_FFOA1)
shapiro.test(aov(M_FFOA1)$residuals)
leveneTest(Average~ Karyotype, FFOA1)


# Histogramm and QQ Plot
model <- aov(M_FFOA1)
par(mar = c(5, 5, 4, 2))
qqPlot(aov(M_FFOA1)$residuals,
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


#Contrast test with reducted data set FFOA1
FFOA1_3RP <- FFOA1[c(which(str_detect(FFOA1$Karyotype, "3RP|STD"))),]
FFOA1_2Lt <- FFOA1[c(which(str_detect(FFOA1$Karyotype, "2Lt|STD"))),]
FFOA1_3RK <- FFOA1[c(which(str_detect(FFOA1$Karyotype, "3RK|STD"))),]

M_FFOA1_3RP <- lm(Average ~ Karyotype, FFOA1_3RP)
averageFFOA1_3RP <- emmeans(M_FFOA1_3RP, "Karyotype")
ComparaisonFFOA1_3RP <- pairs(averageFFOA1_3RP, adjust="none")
ComparaisonFFOA1_3RP <-as.data.frame(ComparaisonFFOA1_3RP)

M_FFOA1_2Lt <- lm(Average ~ Karyotype, FFOA1_2Lt)
averageFFOA1_2Lt <- emmeans(M_FFOA1_2Lt, "Karyotype")
ComparaisonFFOA1_2Lt <- pairs(averageFFOA1_2Lt, adjust="none")
ComparaisonFFOA1_2Lt <-as.data.frame(ComparaisonFFOA1_2Lt)

M_FFOA1_3RK <- lm(Average ~ Karyotype, FFOA1_3RK)
averageFFOA1_3RK <- emmeans(M_FFOA1_3RK, "Karyotype")
ComparaisonFFOA1_3RK <- pairs(averageFFOA1_3RK, adjust="none")
ComparaisonFFOA1_3RK <-as.data.frame(ComparaisonFFOA1_3RK)


Total_contrast <- rbind(ComparaisonFFOA1_3RP,ComparaisonFFOA1_2Lt, ComparaisonFFOA1_3RK )
adjusted_p <- p.adjust(Total_contrast2$p.value, method = "bonferroni")

# alpha adjust
fwa=1-(1-0.05)^(1/18)

#FINAL RESULTS FOR CONTRAST TEST
Significant_Total_contrast <- Total_contrast[Total_contrast$p.value<fwa,] 


# Add p-value
Contrast_FFOA<- separate(Total_contrast, contrast, into= c("group1", "group2"), sep = " - ") 
Contrast_FFOA1 <- mutate(Contrast_FFOA, p.signif=case_when(p.value < fwa ~ "*",TRUE ~ "ns"))   

FFOA_3RP <- FFOA1[c(which(str_detect(FFOA1$Karyotype, "3RP|STD"))),]
FFOA_2Lt <- FFOA1[c(which(str_detect(FFOA1$Karyotype, "2Lt|STD"))),]
FFOA_3RK <- FFOA1[c(which(str_detect(FFOA1$Karyotype, "3RK|STD"))),]
FFOA_HOMO <- FFOA1[-c(which(str_detect(FFOA1$Karyotype, "HET"))),]

p_FFOA_3RP <- Contrast_FFOA1[Contrast_FFOA1$group1 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2") & Contrast_FFOA1$group2 %in% c("STD","3RP","3RP_HET_1","3RP_HET_2"),]
p_FFOA_2Lt <- Contrast_FFOA1[Contrast_FFOA1$group1 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2") & Contrast_FFOA1$group2 %in% c("STD","2Lt","2Lt_HET_1","2Lt_HET_2"),]
p_FFOA_3RK <- Contrast_FFOA1[Contrast_FFOA1$group1 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2") & Contrast_FFOA1$group2 %in% c("STD","3RK","3RK_HET_1","3RK_HET_2"),]
p_FFOA_HOMO <- Contrast_FFOA1[Contrast_FFOA1$group1 %in% c("STD","3RP","2Lt","3RK") & Contrast_FFOA1$group2 %in% c("STD","3RP","2Lt","3RK"),]


# Box Plot
FFOA_B1 <- ggplot(FFOA_3RP, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.525,0.61) +stat_pvalue_manual(p_FFOA_3RP, label="p.signif", hide.ns = TRUE, y.position=0.59, step.increase = 0.08)+
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female femur length 3RP-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

FFOA_B2 <- ggplot(FFOA_3RK, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.525,0.61) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female femur length 3RK-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

FFOA_B3 <- ggplot(FFOA_2Lt, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.525,0.61) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female femur length 2Lt-STD",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))

FFOA_B4 <- ggplot(FFOA_HOMO, aes(x= Karyotype, y = Average, fill = Karyotype)) +
  geom_boxplot() +
  scale_fill_manual(values = Colors, labels = function(x) gsub("_", " ", x)) + scale_x_discrete(labels = function(x) gsub("_", " ", x))+
  ylim(0.525,0.61) +
  theme_bw()+
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.x = element_text(face = "bold", size = 11),
        axis.text.y = element_text(face = "bold", size = 11),
        legend.title = element_text(face = "bold", size = 16),
        legend.text = element_text(face = "bold", size = 12))+
  labs(title = "Female femur length (Homokaryotype)",subtitle = "\u2640", x = "Karyotype", y = expression(bold("Femur length (") * bold(mm) * bold(")"))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(size = 20, face = "bold", hjust = 0))



# General box plot
(FFOA_B1|FFOA_B2)/(FFOA_B3|FFOA_B4) + plot_annotation(title="Female Femur Length", 
                                                      tag_levels='A', 
                                                      theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))



                                                     theme=theme(plot.title=element_text(hjust=0.5, size =18, face="bold"),plot.caption = element_text(hjust=0.5, size =12)))

