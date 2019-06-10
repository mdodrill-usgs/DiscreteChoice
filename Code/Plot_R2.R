###############################################################################
#                                                                     March '19
#        Discrete Choice Stan Model - Rainbow Trout Prey Selection
#                                Plotting
#
#  Notes:
#  * 
#
###############################################################################
setwd("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/")
library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(RColorBrewer)
library(reshape2)

name.key = data.frame(num = c(1:7),
                      name = c("Snails", "Gammarus",
                               "Black Fly Adult", "Black Fly Larva",
                               "Midge Adult", "Midge Larva","Worms"),
                      name2 = c("Snails", "Gammarus",
                                "Black Fly \nAdult / Pupae", "Black Fly \nLarva",
                                "Midge \nAdult / Pupae", "Midge \nLarva","Worms"))

#-----------------------------------------------------------------------------#
# 7 taxa models
# runs = c("Model_Ind_v1_with_R2_Length_2000iter.RData",
#          "Model_Ind_v1_with_R2_Width_2000iter.RData",
#          "Model_Ind_v1_with_R2_Area_2000iter.RData",
#          "Model_Ind_v1_with_R2_Mass_2000iter.RData")

runs = c("Model_Ind_v1_with_R2_Length_2000iter.RData",
         "Model_Ind_v1_with_R2_Width_2000iter.RData",
         "Model_Ind_v1_with_R2_Area_2000iter.RData",
         "Model_Ind_v1_with_R2_Mass_2000iter.RData",
         "Model_Ind_v1_with_R2_Hemi_A_2000iter.RData",
         # "Model_Ind_v1_with_R2_Hemi_p_2000iter.RData",
         "Model_Ind_v1_with_R2_Volume_2000iter.RData")
 
# 6 taxa models
runs = c("Model_Ind_v1_with_R2_Length_2000iter_with_All_RE_no_worm.RData",
         "Model_Ind_v1_with_R2_Width_2000iter_with_All_RE_no_worm.RData",
         "Model_Ind_v1_with_R2_Area_2000iter_with_All_RE_no_worm.RData",
         "Model_Ind_v1_with_R2_Mass_2000iter_with_All_RE_no_worm.RData")

name = "vRE" # this will grab all the v... parameters

new.out = list()

t1 = proc.time()

for(new.i in seq_along(runs)){
  load(runs[new.i])
  
  f1 <- coda::mcmc.list(lapply(1:ncol(fit), function(x) coda::mcmc(as.array(fit)[,x,])))
  f1 <- ggs(f1)
  
  f1$model = rep(runs[new.i], dim(f1)[1])
  f1.sub = f1[which(substr(f1$Parameter,1,nchar(name)) == name),]
  
  f1.sub$R2_re = 1 - (f1.sub[f1.sub$Parameter == "vRE",]$value / f1.sub[f1.sub$Parameter == "vRE_sp_psz_fsz",]$value)
  f1.sub$R2_sp = 1 - (f1.sub[f1.sub$Parameter == "vRE_sp_fsz",]$value / f1.sub[f1.sub$Parameter == "vRE_sp_psz_fsz",]$value)  
  f1.sub$R2_sz = 1 - (f1.sub[f1.sub$Parameter == "vRE_sp",]$value / f1.sub[f1.sub$Parameter == "vRE_sp_psz_fsz",]$value)  
  
  new.out[[new.i]] = f1.sub[,5:8]
  rm(fit)
}

t2 = proc.time()

tmp = do.call(rbind, new.out)

all = melt(tmp, id.vars = "model")

all2 = group_by(all, model, variable) %>%
  summarize(my.mean = mean(value),
            upper = quantile(value, .975),
            lower = quantile(value, .025))

# tmp = rep(c("Area", "Length", "Mass", "Width"), each = 3)

tmp = rep(c("Area", "Hemispherical\n Area", "Length", "Mass",
            "Volume", "Width"), each = 3)

all2$measure = factor(tmp,
                       levels = c("Length", "Area", "Hemispherical\n Area",
                                  "Mass", "Volume", "Width"),
                       ordered = TRUE)

#-----------------------------------------------------------------------------#
# windows(width = 6.5*1.5, height = 5*1.5, record = T, xpos = 25)
windows(width = 3.4*2, height = 3.4*2, record = T, xpos = 25)

all.ltl = all2[all2$variable == "R2_sp",]

# var.lab = ifelse(all.ltl$variable == "R2_sp", "Prey Size", "Prey and Fish Size")
# 
# all.ltl$var.lab = factor(var.lab,
#                          levels = c("Prey Size", "Prey and Fish Size"),
#                          ordered = TRUE) 

p = ggplot(all.ltl, aes(y = my.mean, x = measure)) +
    # geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    # geom_point(position = position_dodge(width = 1), aes(color = measure), size = 3) +
    geom_point(position = position_dodge(width = 1), color = "black", size = 3) +
    # geom_errorbar(aes(ymin = lower, ymax = upper, color = measure),
    geom_errorbar(aes(ymin = lower, ymax = upper), color = "black",
                  position = position_dodge(width = 1),
                  width = 0, size = 1) +
    scale_y_continuous(breaks = c(-.1, 0, .1, .2, .3, .4, .5)) +
    # scale_color_manual(values = c("black", "grey35", "gray50", "gray75")) +
    labs(y = expression(paste('Pseudo R'^' 2',' (+ - 95% CRI)')), x = "", color = "") +
    coord_flip()

p1 = p + theme(axis.title.x = element_text(size = 24, vjust = -.1),
               axis.title.y = element_text(size = 24, vjust = 1),
               axis.text.x = element_text(size = 20, colour = "black"),
               axis.text.y = element_text(size = 20, colour = "black"),
               axis.line = element_line(colour = "black"),
               panel.background = element_rect(fill = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.grid.major = element_line(colour = "white"),
               panel.spacing = unit(.8, "lines"),
               legend.text = element_text(size = 19),
               legend.key = element_rect(fill = "white"),
               legend.position = "none")
               # legend.position = c(.9, .95))
p1

#-----------------------------------------------------------------------------#
# End