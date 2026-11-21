######## statics and figures #####

## 
#libraries 
library(ggplot2)
library(rstatix)
library(dplyr)
library(ggpubr)
library(PMCMRplus)
library(Rmisc)
#loading database 
database_complete <- read.csv("data_base/data_base_complete_YPD_MI_MIF_KCl_0.25_0.5_1_W_NW.csv")
SN_database_complete <- read.csv("data_base/data_base_SN_MI_MIF_KCl_W_NW.csv")
View(database_complete)
View(SN_database_complete)

### sub data base ###
YPD_complete <- database_complete[database_complete$medium == "YPD", 2:5]
View(YPD_complete)
MI_complete <- database_complete[database_complete$medium == "MI", 2:5]
View(MI_complete)
MIF_complete <- database_complete[database_complete$medium == "MIF", 2:5]
View(MIF_complete)
CTRLS <- database_complete[database_complete$condition == "CTRL", 2:5]
View(CTRLS)
SN_CTRLS <- SN_database_complete[SN_database_complete$condition == "CTRL", 2:5]
SN_MIF_OS <- SN_database_complete[SN_database_complete$treatment != "CTRL MI", 2:5]

#### YPD ####
## data visualization ##
g1 <- ggplot(data = YPD_complete, mapping = aes(x = condition, y = area)) +
  geom_boxplot() +
  theme_bw()
g1
g2 <- ggplot(data = YPD_complete, mapping = aes(x = treatment, y = area)) +
  geom_boxplot() +
  theme_bw()
g2
g3 <- ggplot(data = YPD_complete, mapping = aes(x = condition, y = area, 
                                                             colour = treatment)) +
  geom_boxplot() +
  theme_bw()
g3
## descriptive statistics ##
YPD_complete %>% 
  group_by(condition, treatment) %>%
  get_summary_stats(area, type = "mean_sd")

#Normal distribution
YPD_complete %>%
  group_by(condition, treatment) %>%
  shapiro_test(area)

ggqqplot(YPD_complete, "area", ggtheme = theme_bw()) + facet_grid(condition ~ treatment, labeller = "label_both")
#the data has a normal distribution

#homogeneity of variances
bartlett.test(area~ treatment, data = YPD_complete)
bartlett.test(area~ condition, data = YPD_complete)
#There is no homogeneity of variances

#outliers
YPD_complete %>%
  group_by(condition) %>%
  identify_outliers(area)
#there's no outliers

## Statics ##
anova_YPD <- aov(YPD_complete$area ~ YPD_complete$condition * YPD_complete$treatment)
summary(anova_YPD)
#results:
#significant differences in: YPD_complete$condition  p= 0.0355 *
#POST HOC
YPD_complete$condition <- as.factor(YPD_complete$condition)
YPD_complete$treatment <- as.factor(YPD_complete$treatment)
dunnettT3Test(YPD_complete$area ~ YPD_complete$condition)

YPD_combination <- YPD_complete %>%
  mutate(total_group = as.factor(paste(condition, treatment, sep = "_")))
YPD_combination

# post hoc veridic
dunnettT3Test(area ~ total_group, data = YPD_combination)
#Results: inhibition halos of Non washed cells in osmotic shock of 1M KCl
#are significant different of inhibition halos of control cells
#p = 0.014

### Figure 1b ### # # #
#Needed colors palette
Col_YPD <- c("#E08736", "thistle", "#C53C8F")
col2rgb(Col_YPD)
#colores_YPD <- t(col2rgb(colores_YPD)) #to convert and maintain as 8-bit color as required
data_graph_YPD <- summarySE(YPD_complete, measurevar = "area", groupvars = c("condition", "treatment"))
data_graph_YPD
data_graph_YPD$condition <- factor(data_graph_YPD$condition, 
                                 levels = c("CTRL", "W", "NW"))
graph_Fig1b_YPD <- ggplot(data_graph_YPD, aes(x = treatment, y = area, fill = condition))+
  geom_bar(position = position_dodge(0.96), stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = Col_YPD, name = "Condition", labels = c("Untreated", "Washed cells", 
                                                                           "Non\nwashed cells"))
graph_Fig1b_YPD  
graph_Fig1b_YPD2 <- graph_Fig1b_YPD + ggtitle("Inhibition halos produced after\nosmotic shock",
                                              subtitle = "Cells grown in YPD medium")+
  
  labs(x = "Treatment", y = bquote('Area'~(cm^2)))+
  theme(plot.title = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2.5)),
        plot.subtitle = element_text(hjust = 0.5, lineheight = .8, size = rel(1.7)))+
  theme(axis.title.x = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.title.y = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.text.x = element_text(vjust = 0.5, size = rel(1)))+
  theme(panel.background = element_rect(fill = 'gray98'), panel.grid = element_line(color = "gray90"))+
  theme(axis.text = element_text(size = 17, face = "bold", vjust = 0.5))+
  theme(legend.text = element_text(size = 12, face = "bold"))+
  theme(legend.title = element_text(size = 14, face = "bold"))+
  scale_x_discrete(labels = c("CTRL YPD" = expression(bold("Untreated")), "KCl 0.25" = expression(bold("KCl 0.25 M")),
                              "KCl 0.5" = expression(bold("KCl 0.5 M")), "KCl 1" = expression(bold("KCl 1 M"))),
                   guide = guide_axis(angle = 25))+
  geom_signif(y_position = 0.43, xmin = 1, xmax = 4.23, annotations = "*",
              tip_length = 0.05, textsize=8, vjust = 0.4)
graph_Fig1b_YPD2

#### MI, MIF, YPD CTRLS ####
## data visualization ##
attach(CTRLS)
g1 <- ggplot(data = CTRLS, mapping = aes(x = treatment, y = area)) +
  geom_boxplot() +
  theme_bw()
g1

## descriptive statistics ##
CTRLS %>% 
  group_by(treatment) %>%
  get_summary_stats(area, type = "mean_sd")

#Normal distribution
CTRLS %>%
  group_by(treatment) %>%
  shapiro_test(area)
ggqqplot(CTRLS, "area", ggtheme = theme_bw()) + facet_grid(medium ~ treatment, labeller = "label_both")
#the data has a normal distribution

#homogeneity of variances
bartlett.test(area~ treatment, data = CTRLS)
bartlett.test(area~ medium, data = CTRLS)
#there's no homogeneity of variances

#outliers
CTRLS %>%
  group_by(treatment) %>%
  identify_outliers(area)
#there's no outliers

## Statics ##
anova_CTRLS <- aov(CTRLS$area ~ CTRLS$treatment)
summary(anova_CTRLS)
#results:
#significant differences in: CTRLS$treatment  p= 3.57e-12 ***
#POST HOC
CTRLS$treatment <-  as.factor(CTRLS$treatment)
dunnettT3Test(CTRLS$area ~ CTRLS$treatment)
#Results: inhibition halos produced by cells grown in YPD medium are significant different than those
#inhibition halos produced by cells grown in MI and MIF media, but not between MI and MIF halos
#YPD - MI p =1.4e-12
#YPD - MIF p = 4.3e-08 
#MI - MIF p= 0.96

### Figure 1c ### # # #
#Needed colors palette
Col_YPD_MI <- c("#E08736", "goldenrod1")
col2rgb(Col_YPD_MI)
YPD_MI_CTRL <- CTRLS[1:32, 1:4]
YPD_MI_CTRL
data_graph_YPD_MI <- summarySE(YPD_MI_CTRL, measurevar = "area", groupvars = c("treatment"))
data_graph_YPD_MI
data_graph_YPD_MI$treatment <- factor(data_graph_YPD_MI$treatment, 
                                   levels = c("CTRL YPD", "CTRL MI"))

graph_Fig1c_YPD_MI <- ggplot(data_graph_YPD_MI, aes(x = treatment, y = area, fill = treatment))+
  geom_bar(position = position_dodge(0.96), stat = "identity", colour = "black", width = 0.65)+
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = Col_YPD_MI)
graph_Fig1c_YPD_MI 

graph_Fig1c_YPD_MI2 <- graph_Fig1b_YPD_MI + ggtitle("Inhibition halos produced by cells\ngrown YPD and MI media")+
  labs(x = "Medium", y = bquote('Area'~(cm^2)))+
  theme(plot.title = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2.5)))+
  theme(axis.title.x = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.title.y = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.text.x = element_text(vjust = 0.5, size = rel(1)))+
  theme(panel.background = element_rect(fill = 'gray98'), panel.grid = element_line(color = "gray90"))+
  theme(axis.text = element_text(size = 17, face = "bold", vjust = 0.5))+
  geom_signif(y_position = 0.5, xmin = 1, xmax = 2, annotations = "***",
              tip_length = 0.05, textsize=8, vjust = 0.4)+
  guides(fill="none")+
  scale_x_discrete(labels = c("CTRL YPD" = expression(bold("YPD")), "CTRL MI" = expression(bold("MI"))))+
  theme(plot.tag = element_text(face = "bold", size = 20))
graph_Fig1c_YPD_MI2

### Figure 2b ### # # #
#Needed colors palette
Col_CTRLS <- c("#E08736", "goldenrod1", "tomato3")
col2rgb(Col_CTRLS)
data_graph_CTRLS <- summarySE(CTRLS, measurevar = "area", groupvars = c("treatment"))
data_graph_CTRLS

data_graph_CTRLS$treatment <- factor(data_graph_CTRLS$treatment, 
                                      levels = c("CTRL YPD", "CTRL MI", "CTRL MIF"))

graph_Fig2b_CTRLS <- ggplot(data_graph_CTRLS, aes(x = treatment, y = area, fill = treatment))+
  geom_bar(position = position_dodge(0.96), stat = "identity", colour = "black", width = 0.75)+
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = Col_CTRLS)
graph_Fig2b_CTRLS
graph_Fig2b_CTRLS2 <- graph_Fig2b_CTRLS + ggtitle("Inhibition halos produced by\ncells grown in different media")+
  labs(x = "Medium", y = bquote('Area'~(cm^2)))+
  theme(plot.title = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2.3)))+
  theme(axis.title.x = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.title.y = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.text.x = element_text(vjust = 0.5, size = rel(1)))+
  theme(panel.background = element_rect(fill = 'gray98'), panel.grid = element_line(color = "gray90"))+
  theme(axis.text = element_text(size = 17, face = "bold", vjust = 0.5))+
  geom_signif(y_position = c(0.5, 0.6), xmin = c(1, 1), xmax = c(2, 3), annotations = c("***", "***"),
              tip_length = 0.05, textsize=8, vjust = 0.4)+
  geom_signif(y_position = 0.55, xmin = 2, xmax = 3, annotations = "NS", tip_length = 0.05, textsize = 5, 
              vjust = 0.2)+
  guides(fill="none")+
  scale_x_discrete(labels = c("CTRL YPD" = expression(bold("YPD")), "CTRL MI" = expression(bold("MI")), 
                              "CTRL MIF" = expression(bold("MIF"))))
  
graph_Fig2b_CTRLS2


#### Osmotic shock MI ####
detach(CTRLS)
attach(MI_complete)
## data visualization ##
g1 <- ggplot(data = MI_complete, mapping = aes(x = condition, y = area)) +
  geom_boxplot() +
  theme_bw()
g1
g2 <- ggplot(data = MI_complete, mapping = aes(x = treatment, y = area)) +
  geom_boxplot() +
  theme_bw()
g2
g3 <- ggplot(data = MI_complete, mapping = aes(x = condition, y = area, 
                                                colour = treatment)) +
  geom_boxplot() +
  theme_bw()
g3

## descriptive statistics ##
MI_complete %>% 
  group_by(condition, treatment) %>%
  get_summary_stats(area, type = "mean_sd")

#Normal distribution
MI_complete %>%
  group_by(condition, treatment) %>%
  shapiro_test(area)
ggqqplot(MI_complete, "area", ggtheme = theme_bw()) + facet_grid(condition ~ treatment, labeller = "label_both")
#the data has a normal distribution

#homogeneity of variances
bartlett.test(area~ treatment, data = MI_complete)
bartlett.test(area~ condition, data = MI_complete)
#There is no homogeneity of variances

#outliers
MI_complete %>%
  group_by(condition, treatment) %>%
  identify_outliers(area)
#there's no outliers

## Statics ##
anova_MI_OS <- aov(MI_complete$area ~  MI_complete$condition * MI_complete$treatment)
summary(anova_MI_OS)
#results:
#significant differences in: MI_complete$treatment  p= 2.77e-05 ***
#no significant differences in: MI_complete$condition  p=   0.326 
#POST HOC
MI_complete$treatment <- as.factor(MI_complete$treatment)
MI_complete$condition <- as.factor(MI_complete$condition)
#all the posible combinations
MI_complete_combination <- MI_complete %>%
  mutate(total_group = as.factor(paste(condition, treatment, sep = "_")))
MI_complete_combination
# post hoc veridic
dunnettT3Test(area ~ total_group, data = MI_complete_combination)
#results:inhibition halos of Non washed and washed cells in osmotic shock of 1M KCl
#are significant different of inhibition halos of control cells (MI medium):
#CTRL - W_KCl p= 0.00106
#CTRL - NW_KCl p= 0.00017
#but not between them
#NW_KCl - W_KCl p =  0.78107

### Figure 1d ### # # #
#Needed colors palette
Col_MI <- c("goldenrod1", "dodgerblue4", "deepskyblue3")
col2rgb(Col_MI)

data_graph_MI_OS <- summarySE(MI_complete, measurevar = "area", groupvars = c("condition", "treatment"))
data_graph_MI_OS
data_graph_MI_OS$condition <- factor(data_graph_MI_OS$condition, 
                                   levels = c("CTRL", "W", "NW"))
graph_Fig1d_MI_OS <- ggplot(data_graph_MI_OS, aes(x = treatment, y = area, fill = condition))+
  geom_bar(position = position_dodge(0.96), stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = Col_MI, name = "Condition", labels = c("Untreated", "Washed cells", 
                                                                     "Non\nwashed cells"))
graph_Fig1d_MI_OS  
graph_Fig1d_MI_OS2 <- graph_Fig1d_MI_OS + ggtitle("Inhibition halos produced after\nosmotic shock",
                                                  subtitle = "Cells grown in MI medium")+
  labs(x = "Treatment", y = bquote('Area'~(cm^2)))+
  theme(plot.title = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2.5)),
        plot.subtitle = element_text(hjust = 0.5, lineheight = .8, size = rel(1.7)))+
  theme(axis.title.x = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.title.y = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.text.x = element_text(vjust = 0.5, size = rel(1)))+
  theme(panel.background = element_rect(fill = 'gray98'), panel.grid = element_line(color = "gray90"))+
  theme(axis.text = element_text(size = 17, face = "bold", vjust = 0.5))+
  geom_signif(y_position = 0.8, xmin = 1, xmax = 2.225, annotations = "***",
              tip_length = 0.05, textsize=8, vjust = 0.4)+
  geom_signif(y_position = 0.72, xmin = 1.775, xmax = 2.225, annotations = "NS",  tip_length = 0.05, textsize = 5,
              vjust = 0.2)+
  theme(legend.text = element_text(size = 12, face = "bold"))+
  theme(legend.title = element_text(size = 14, face = "bold"))+
  scale_x_discrete(labels = c("CTRL MI" = expression(bold("Untreated")), "KCl" = expression(bold("KCl"))))+
  theme(plot.tag = element_text(face = "bold", size = 20))
graph_Fig1d_MI_OS2



#### Inhibition halos by supernatants ####
## data visualization ##
detach(MI_complete)
attach(SN_CTRLS)
g1 <- ggplot(data = SN_CTRLS, mapping = aes(x = treatment, y = area)) +
  geom_boxplot() +
  theme_bw()
g1

## descriptive statistics ##
SN_CTRLS %>% 
  group_by(treatment) %>%
  get_summary_stats(area, type = "mean_sd")

#Normal distribution
SN_CTRLS %>%
  group_by(treatment) %>%
  shapiro_test(area)
ggqqplot(SN_CTRLS, "area", ggtheme = theme_bw()) + facet_grid(medium ~ treatment, labeller = "label_both")
#the data has a normal distribution

#homogeneity of variances
bartlett.test(area~ treatment, data = SN_CTRLS)
bartlett.test(area~ medium, data = SN_CTRLS)
#There's no homogeneity of variances

#outliers
SN_CTRLS %>%
  group_by(treatment) %>%
  identify_outliers(area)
#there's one outlier, no extreme.

## Statics ##
t.test(SN_CTRLS[SN_CTRLS$treatment == "CTRL MI", 4],
       SN_CTRLS[SN_CTRLS$treatment == "CTRL MIF", 4], paired = FALSE, var.equal = FALSE)
#Result:
#significant differences in between MI and MIF inhibitions halos produced by supernatant
# p-value < 2.2e-16

### Figure 2c ### # # #
#Needed colors palette
Col_SN_CTRLS <- c( "goldenrod1", "tomato3")

data_graph_SN_CTRLS <- summarySE(SN_CTRLS, measurevar = "area", groupvars = c("treatment"))
data_graph_SN_CTRLS

graph_Fig2c_SN_CTRLS <- ggplot(data_graph_SN_CTRLS, aes(x = treatment, y = area, fill = treatment))+
  geom_bar(position = position_dodge(0.96), stat = "identity", colour = "black", alpha = 0.5, width = 0.6)+
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = Col_SN_CTRLS)
graph_Fig2c_SN_CTRLS
graph_Fig2c_SN_CTRLS2 <- graph_Fig2c_SN_CTRLS + ggtitle("Inhibition halos produced by supernatants of\ncells grown in different media")+
  labs(x = "Medium", y = bquote('Area'~(cm^2)))+
  theme(plot.title = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2.3)))+
  theme(axis.title.x = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.title.y = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.text.x = element_text(vjust = 0.5, size = rel(1)))+
  theme(panel.background = element_rect(fill = 'gray98'), panel.grid = element_line(color = "gray90"))+
  theme(axis.text = element_text(size = 17, face = "bold", vjust = 0.5))+
  geom_signif(y_position = 0.67, xmin = 1, xmax = 2, annotations = "***",
              tip_length = 0.03, textsize=8, vjust = 0.4)+
  guides(fill="none")+
  scale_x_discrete(labels = c("CTRL MI" = expression(bold("MI")), "CTRL MIF" = expression(bold("MIF"))))+
  theme(plot.tag = element_text(face = "bold", size = 20))
graph_Fig2c_SN_CTRLS2


#### Osmotic shock MIF ####
detach(SN_CTRLS)
attach(MIF_complete)

## data visualization ##
g1 <- ggplot(data = MIF_complete, mapping = aes(x = condition, y = area)) +
  geom_boxplot() +
  theme_bw()
g1
g2 <- ggplot(data = MIF_complete, mapping = aes(x = treatment, y = area)) +
  geom_boxplot() +
  theme_bw()
g2
g3 <- ggplot(data = MIF_complete, mapping = aes(x = condition, y = area, 
                                               colour = treatment)) +
  geom_boxplot() +
  theme_bw()
g3

## descriptive statistics ##
MIF_complete %>% 
  group_by(condition, treatment) %>%
  get_summary_stats(area, type = "mean_sd")

#Normal distribution
MIF_complete %>%
  group_by(condition, treatment) %>%
  shapiro_test(area)
#the data has a normal distribution

#homogeneity of variances
bartlett.test(area~ treatment, data = MIF_complete)
bartlett.test(area~ condition, data = MIF_complete)
#There IS homogeneity of variances

#outliers
MIF_complete %>%
  group_by(condition, treatment) %>%
  identify_outliers(area)
#there's one outlier, no extreme.

## Statics ##
anova_MIF_OS <- aov(MIF_complete$area ~ MIF_complete$condition * MIF_complete$treatment)
summary(anova_MIF_OS)
#POST HOC
TukeyHSD(anova_MIF_OS)
#results:
#significant differences in: W-CTRL  p = 0.0000621
#significant differences in: W-NW    p = 0.0280973
# no significant differences in: NW-CTRL p = 0.1050708

### Figure 2d ### # # #
#Needed colors palette
Col_MIF <- c("tomato3", "palegreen4", "darkolivegreen3")
col2rgb(Col_MIF)

data_graph_MIF_OS <- summarySE(MIF_complete, measurevar = "area", groupvars = c("condition", "treatment"))
data_graph_MIF_OS
data_graph_MIF_OS$condition <- factor(data_graph_MIF_OS$condition, 
                                     levels = c("CTRL", "W", "NW"))
graph_Fig2d_MIF_OS <- ggplot(data_graph_MIF_OS, aes(x = treatment, y = area, fill = condition))+
  geom_bar(position = position_dodge(0.96), stat = "identity", colour = "black")+
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = Col_MIF, name = "Condition", labels = c("Untreated", "Washed cells", 
                                                                    "Non\nwashed cells"))
graph_Fig2d_MIF_OS  
graph_Fig2d_MIF_OS2 <- graph_Fig2d_MIF_OS + ggtitle("Inhibition halos produced after\nosmotic shock",
                                                    subtitle = "Cells grown in MIF medium")+
  labs(x = "Treatment", y = bquote('Area'~(cm^2)))+
  theme(plot.title = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2.3)),
        plot.subtitle = element_text(hjust = 0.5, lineheight = .8, size = rel(1.7)))+
  theme(axis.title.x = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.title.y = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.text.x = element_text(vjust = 0.5, size = rel(1)))+
  theme(panel.background = element_rect(fill = 'gray98'), panel.grid = element_line(color = "gray90"))+
  theme(axis.text = element_text(size = 17, face = "bold", vjust = 0.5))+
  geom_signif(y_position = c(0.78, 0.73), xmin = c(1, 1.775), xmax = c(1.775, 2.225), annotations = c("***", "*"),
              tip_length = 0.05, textsize=8, vjust = 0.4)+
  theme(legend.text = element_text(size = 12, face = "bold"))+
  theme(legend.title = element_text(size = 14, face = "bold"))+
  scale_x_discrete(labels = c("CTRL MIF" = expression(bold("Untreated")), "KCl" = expression(bold("KCl"))))+
  theme(plot.tag = element_text(face = "bold", size = 20))
graph_Fig2d_MIF_OS2


#### Inhibtion halos by supernantans in Osmotic shock MIF ####
detach(MIF_complete)
attach(SN_MIF_OS)

g1 <- ggplot(data = SN_MIF_OS, mapping = aes(x = condition, y = area)) +
  geom_boxplot() +
  theme_bw()
g1
g2 <- ggplot(data = SN_MIF_OS, mapping = aes(x = treatment, y = area)) +
  geom_boxplot() +
  theme_bw()
g2
g3 <- ggplot(data = SN_MIF_OS, mapping = aes(x = condition, y = area, 
                                                colour = treatment)) +
  geom_boxplot() +
  theme_bw()
g3

## descriptive statistics ##
SN_MIF_OS %>% 
  group_by(condition, treatment) %>%
  get_summary_stats(area, type = "mean_sd")

#Normal distribution
SN_MIF_OS %>%
  group_by(condition, treatment) %>%
  shapiro_test(area)
#the data has a normal distribution

#homogeneity of variances
bartlett.test(area~ treatment, data = SN_MIF_OS)
bartlett.test(area~ condition, data = SN_MIF_OS)
#There's no homogeneity of variances

#outliers
SN_MIF_OS %>%
  group_by(condition, treatment) %>%
  identify_outliers(area)
#there's four outliers, no extreme.

## Statics ##
anova_SN_MIF_OS <- aov(SN_MIF_OS$area ~ SN_MIF_OS$condition * SN_MIF_OS$treatment)
summary(anova_SN_MIF_OS)
#significant differences in: SN_MIF_OS$condition  p <2e-16 ***
SN_MIF_OS$treatment <- as.factor(SN_MIF_OS$treatment)
SN_MIF_OS$condition <- as.factor(SN_MIF_OS$condition)
SN_MIF_complete_combination <- SN_MIF_OS %>%
  mutate(total_group = as.factor(paste(condition, treatment, sep = "_")))
SN_MIF_complete_combination
# post hoc veridic
dunnettT3Test(area ~ total_group, data = SN_MIF_complete_combination)
#significant differences in: 
#CTRL - NW_KCl p = 3.7e-08
#CTRL - W_KCl  p < 2e-16
#NW_KCl - W_KCl  p = 0.019

### Figure 2e ### # # #
data_graph_SN_MIF_OS <- summarySE(SN_MIF_OS, measurevar = "area", groupvars = c("condition", "treatment"))
data_graph_SN_MIF_OS
data_graph_SN_MIF_OS$condition <- factor(data_graph_MIF_OS$condition, 
                                      levels = c("CTRL", "W", "NW"))
graph_Fig2e_SN_MIF_OS <- ggplot(data_graph_SN_MIF_OS, aes(x = treatment, y = area, fill = condition))+
  geom_bar(position = position_dodge(0.96), stat = "identity", colour = "black", alpha = 0.5)+
  geom_errorbar(aes(ymin=area-se, ymax=area+se), width=0.2, position = position_dodge(0.9))+
  scale_fill_manual(values = Col_MIF, name = "Condition", labels = c("Untreated", "Washed cells", 
                                                                     "Non\nwashed cells"))
graph_Fig2e_SN_MIF_OS 
graph_Fig2e_SN_MIF_OS2 <- graph_Fig2e_SN_MIF_OS + ggtitle("Inhibition halos produced by supernatant\nafter osmotic shock",
                                                          subtitle = "Cells grown in MIF medium")+
  labs(x = "Treatment", y = bquote('Area'~(cm^2)))+
  theme(plot.title = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2.3)), 
        plot.subtitle = element_text(hjust = 0.5, lineheight = .8, size = rel(1.7)))+
  theme(axis.title.x = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.title.y = element_text(hjust = 0.5, lineheight = .8, face = "bold", size = rel(2)))+
  theme(axis.text.x = element_text(vjust = 0.5, size = rel(1)))+
  theme(panel.background = element_rect(fill = 'gray98'), panel.grid = element_line(color = "gray90"))+
  theme(axis.text = element_text(size = 17, face = "bold", vjust = 0.5))+
  geom_signif(y_position = c(0.24, 0.22, 0.26), xmin = c(1, 1.775, 1), xmax = c(1.775, 2.225, 2.225), annotations = c("***", "*", "***"),
              tip_length = 0.025, textsize=8, vjust = 0.4)+
  theme(legend.text = element_text(size = 12, face = "bold"))+
  theme(legend.title = element_text(size = 14, face = "bold"))+
  scale_x_discrete(labels = c("CTRL MIF" = expression(bold("Untreated")), "KCl" = expression(bold("KCl"))))+
  theme(plot.tag = element_text(face = "bold", size = 20))
graph_Fig2e_SN_MIF_OS2

### Every figure was manually saved in personal PC "Figures carpet" ###

## For  complete figures including photographs of halos #####
# ##Fig1 ## ## ##
install.packages("cowplot")
library(cowplot)
Figure_1a <- ggdraw() + 
  draw_image("Figures/Fig1a_modified.png", scale = 0.95)
Figure_1a_new <- ggdraw() +
  draw_image("Figures/Fig1a_new_modified.png", scale = 0.95)

tema_paneles <- theme(
  plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
  plot.subtitle = element_text(size = 15, hjust = 0.5),
  plot.title.position = "plot", # Centra el título en todo el ancho del panel
  plot.margin = margin(t = 22, r = 10, b = 10, l = 10) # Espacio arriba para la letra a, b, c...
)
graph_Fig1b_YPD2    <- graph_Fig1b_YPD2 + tema_paneles +  
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.10)))
graph_Fig1c_YPD_MI2 <- graph_Fig1c_YPD_MI2 + tema_paneles + 
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.18)))
graph_Fig1d_MI_OS2  <- graph_Fig1d_MI_OS2 + tema_paneles + 
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.18)))

#uniendo paneles superiores 
fila_superior <- plot_grid(
  Figure_1a, Figure_1a_new, graph_Fig1b_YPD2,
  labels = c("a", "b", "c"),
  label_size = 18,
  label_fontface = "bold",
  ncol = 3,
  rel_widths = c(1, 1, 1.15) # Da un poco más de ancho al panel c por su leyenda
)
#uniendo paneles inferiores
fila_inferior <- plot_grid(
  NULL, graph_Fig1c_YPD_MI2, graph_Fig1d_MI_OS2, NULL,
  labels = c("", "d", "e", ""),
  label_size = 18,
  label_fontface = "bold",
  ncol = 4,
  rel_widths = c(0.4, 1, 1.15, 0.4) # Los 0.4 a los lados centran ambos gráficos
)

#todos los pabeles en la figura 1
Figure1_final <- plot_grid(
  fila_superior,
  fila_inferior,
  ncol = 1,
  rel_heights = c(1, 1)
)

Figure1_final
#guardando figura 1 completa 

ggsave(
  filename = "Figures/figure1_complete_modified.pdf",
  plot = Figure1_final,
  width = 16,
  height = 9,
  units = "in"
)

ggsave(
  filename = "Figures/figure1_complete_modified.jpg",
  plot = Figure1_final,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300
)

# ## Fig2 ## ## ##
Figure_2a <- ggdraw() + 
  draw_image("Figures/Fig2a_modified.png", scale = 0.95)

graph_Fig2b_CTRLS2 <- graph_Fig2b_CTRLS2 + tema_paneles +  
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.10)))
graph_Fig2c_SN_CTRLS2 <-graph_Fig2c_SN_CTRLS2 + tema_paneles + 
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.18)))
graph_Fig2d_MIF_OS2 <- graph_Fig2d_MIF_OS2 + tema_paneles + 
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.18)))
graph_Fig2e_SN_MIF_OS2 <- graph_Fig2e_SN_MIF_OS2 + tema_paneles + 
  scale_y_continuous(expand = expansion(mult = c(0.09, 0.18)))
#superior
Figure2_superior <- plot_grid(
  Figure_2a, graph_Fig2b_CTRLS2, graph_Fig2c_SN_CTRLS2,
  labels = c("a", "b", "c"),
  label_size = 18,         
  label_fontface = "bold",
  ncol = 3,
  rel_widths = c(1, 1, 1.15)  
)
Figure2_superior
#inferior
Figure2_inferior <- plot_grid(
  NULL, graph_Fig2d_MIF_OS2, graph_Fig2e_SN_MIF_OS2, NULL,
  labels = c("", "d", "e", ""),
  label_size = 18,         
  label_fontface = "bold",
  ncol = 4,
  rel_widths = c(0.4, 1, 1.15, 0.4)
)
Figure2_inferior
#todos los paneles juntos
Figure2_final <- plot_grid(
  Figure2_superior,
  Figure2_inferior,
  ncol = 1,
  rel_heights = c(1, 1)
)
Figure2_final

#guardando figura completa ##
ggsave(
  filename = "Figures/figure2_complete_modified.pdf",
  plot = Figure2_final,
  width = 16,
  height = 9,
  units = "in"
)

ggsave(
  filename = "Figures/figure2_complete_modified.jpg",
  plot = Figure2_final,
  width = 16,
  height = 9,
  units = "in",
  dpi = 300
)





ggsave(
  filename = "Figures/figure2_complete.pdf",
  plot = Figure2_final,
  width = 15,
  height = 9,
  units = "in"
)

ggsave(
  filename = "Figures/figure2_complete.jpg",
  plot = Figure2_final,
  width = 15,
  height = 9,
  units = "in",
  dpi = 300
)
