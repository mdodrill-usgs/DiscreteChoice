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
library(fishR)
library(gridExtra)
library(ggthemes)
library(ggplot2)
library(dplyr)
library(plot3D)
library(tidyr)

load("Model_Ind_v1_with_R2_Width_2000iter.RData")
fit.1 = fit
rm(fit)

name.key = data.frame(num = c(1:7),
                      name = c("Snails", "Gammarus",
                               "Black Fly Adult", "Black Fly Larva",
                               "Midge Adult", "Midge Larva","Worms"),
                      name2 = c("Snails", "Gammarus",
                                "Black Fly \nAdults", "Black Fly \nLarvae",
                                "Midge \nAdults", "Midge \nLarvae","Worms"),
                      name3 = c("Snails", "Gammarus",
                                "Black Fly Adult / Pupae", "Black Fly Larva",
                                "Midge Adult / Pupae", "Midge Larva","Worms"))

#-----------------------------------------------------------------------------#
# these are the 'fishR' tools 
# just plotting the terms...

runs = c("Model_Ind_v1_with_R2_Length_2000iter.Rdata",
         "Model_Ind_v1_with_R2_Width_2000iter.Rdata",
         "Model_Ind_v1_with_R2_Area_2000iter.Rdata",
         "Model_Ind_v1_with_R2_Mass_2000iter.Rdata")

names = c("Length", "Width", "Area", "Mass")

outter = list()

for(ii in 1:4){
  load(runs[ii])
  
  tmp = organize(fit, par.name = "beta_f_sz_int")
  
  all = group_by(tmp, Parameter) %>%
    summarize(my.mean = mean(exp(value)),
              upper = quantile(exp(value), .975),
              lower = quantile(exp(value), .025))
  
  all$name = names[ii]
  
  outter[[ii]] = all
  rm(fit)
}

big = do.call(rbind, outter)

big$name = factor(big$name, 
                  levels = c("Length", "Width", "Area", "Mass"),
                  ordered = TRUE)

#--------------------------------------
windows(width = 6.5, height = 5*1.5, record = T, xpos = 25)

p = ggplot(big, aes(x = name, y = my.mean)) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, size = 2) +
  labs(x = "", y = "Fish Size Effect (+ - 95% CRI)") +
  theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  labs(title = "C")

g = p + theme(plot.margin = unit(c(1,1,1,1), "cm"),
              panel.spacing = unit(1, "lines"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title.y = element_text(size = 28, vjust = 2),
              axis.ticks = element_line(color = "black"),
              title = element_text(size = 24),
              axis.text.y = element_text(size = 24, colour = "black"),
              axis.text.x = element_text(size = 22, colour = "black"))#,
g

#-----------------------------------------------------------------------------#
# plotting the mean fish size relationship with some uncertanty...

d = bayes_summary(fit.1, par.name = "beta_f_sz_int")
a = bayes_summary(fit.1, par.name = "mu_sp_all")
b = bayes_summary(fit.1, par.name = "beta_sz")  

a.1 = as.numeric(a[6,2])
b.1 = as.numeric(b[1,2])
d.1 = as.numeric(d[1,2])

size = seq(min(fish_sz), max(fish_sz), 1)
#--------------------------------------
tmp.1 = organize(fit = fit.1, par.name = "beta_f_sz_int")
tmp.2 = organize(fit = fit.1, par.name = "beta_sz")
tmp.2 = tmp.2[tmp.2$Parameter == "beta_sz",]

tmp.a = organize(fit.1, par.name = "mu_sp_all")
tmp.a$id = as.numeric(gsub("[^0-9]", "", as.character(tmp.a$Parameter)))
tmp.a.sp = tmp.a[tmp.a$id == 6,]


a.1 = as.numeric(a[6,2])
b.1 = as.numeric(b[1,2])

out.list = list()

for(i in 1:500){
  # tmp.d = a.1 + b.1 * as.numeric(tmp.1[i,4]) * size # uncertainty in the fish slope part
  tmp.d = as.numeric(tmp.a.sp[i,4]) + b.1 * as.numeric(tmp.1[i,4]) * size # uncertainty in the fish slope part, and prey int
  # tmp.d = a.1 + as.numeric(tmp.2[i,4]) * as.numeric(tmp.1[i,4]) * size
  
  out.list[[i]] = data.frame(run = i,
                             fish_len = size,
                             fish_lab = seq(min(t.sz, na.rm = T), max(t.sz, na.rm = T),1),
                             pred_log = tmp.d,
                             pred = exp(tmp.d))
}


all.out = do.call(rbind, out.list)

mean.tmp = tmp = a.1 + b.1 * d.1 * size
mean.d = data.frame(run = 1,
                    fish_len = size,
                    fish_lab = seq(min(t.sz, na.rm = T), max(t.sz, na.rm = T),1),
                    pred_log = mean.tmp,
                    pred = exp(mean.tmp))


#--------------------------------------
p3 = ggplot(all.out, aes(x = fish_lab, y = pred, group = run)) +
  geom_line(alpha = .2, color = "gray") +
  geom_line(data = mean.d, aes(x = fish_lab, y = pred, group = run),
            color = "black", size = 2) +
  labs(x = "Fork Length (mm)", y = "Selection for Reference \n Prey Type and Size",
       title = "") +
  theme_minimal() +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) +
  coord_cartesian(ylim = c(0,300)) +
  labs(title = "B")

g = p3 + theme(plot.margin = unit(c(1,1,1,1), "cm"),
               panel.spacing = unit(1, "lines"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               axis.title.y = element_text(size = 28, vjust = 4),
               axis.title.x = element_text(size = 28, vjust = -2),
               title = element_text(size = 24),
               axis.ticks = element_line(color = "black"),
               axis.text.y = element_text(size = 24, colour = "black"),
               axis.text.x = element_text(size = 28, colour = "black"))#,
g


#-----------------------------------------------------------------------------#
# 3D Plot! 
windows(width = 6.5*1.5, height = 5*1.5, record = T, xpos = 25)

d = bayes_summary(fit.1, par.name = "beta_f_sz_int")
a = bayes_summary(fit.1, par.name = "mu_sp_all")
b = bayes_summary(fit.1, par.name = "beta_sz")  

a.1 = as.numeric(a[6,2])  # midge larvae
b.1 = as.numeric(b[1,2])
d.1 = as.numeric(d[1,2])

size = seq(min(fish_sz), max(fish_sz), 1)

p.sz = sz[spp == 6]  # midge larvae

# cut this down a little for plotting
p.sz = seq(min(p.sz), max(p.sz), .1)[40:60]

out = list()

for(i in 1:length(p.sz)){
  tmp = a.1 + b.1 * p.sz[i] + d.1 * size  
  
  out[[i]] =  data.frame(fish_len = size,
                         fish_lab = seq(min(t.sz, na.rm = T), max(t.sz, na.rm = T),1),
                         pred_log = tmp,
                         pred = exp(tmp),
                         p.sz = p.sz[i])
}

all.out = do.call(rbind, out)

# format for presp3D
test = select(all.out, fish_lab, p.sz, pred) %>%
       spread(p.sz, pred)

test2 = as.matrix(test[,-1])

# standardize
test3 = test2 / mean(test2)

# cut down a little for plotting
sm = test3[seq(1, nrow(test3), by = 10),
           seq(1, ncol(test3), by = 1)]

persp3D(z = sm, col = ramp.col(c("white", "blue")), border = "black",
        bty = "b2", ticktype = "detailed", #resfac = 2,
        xlab = "Fish Size", ylab = "Prey Size",
        zlab = "Relative Selection",
        theta = 60, phi = 2, clab = NULL, cex.lab = 2)

# axis(1, at = seq(0, 1, by = .1))

#-----------------------------------------------------------------------------#
# library(rgl)
# plot3d(x = all.out$fish_lab, y = all.out$pred, z = all.out$p.sz)
#-----------------------------------------------------------------------------#
# make the plot with custom labels (see working_3D.R)

par(xpd = TRUE)

pmat = persp3D(z = sm, col = ramp.col(c("white", "blue")), border = "black",
               bty = "b2",# ticktype = "detailed", #resfac = 2,
               xlab = "Fish Size", ylab = "Prey Size",
               zlab = "Relative Selection",
               theta = 60, phi = 2, clab = NULL, cex.lab = 2,
               axes = FALSE)

x.axis <- seq(0,1,.2)
min.x <- 0
max.x <- 1
y.axis <- seq(0,1,.25)
min.y <- 0
max.y <- 1
z.axis <- seq(0, 10, by=2)
min.z <- 0
max.z <- 10

tick.start <- trans3d(x.axis, min.y, min.z, pmat)
tick.end <- trans3d(x.axis, (min.y - 0.05), min.z, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

tick.start <- trans3d(max.x, y.axis, min.z, pmat)
tick.end <- trans3d(max.x + 0.05, y.axis, min.z, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)

tick.start <- trans3d(min.x, min.y, z.axis, pmat)
tick.end <- trans3d(min.x, (min.y - 0.05), z.axis, pmat)
segments(tick.start$x, tick.start$y, tick.end$x, tick.end$y)


labels <- c("50", "120", "190", "260", "330", "400")
label.pos <- trans3d(x.axis, (min.y - 0.15), min.z, pmat)
text(label.pos$x, label.pos$y, labels = labels, adj=c(0, NA), srt=0, cex=1.5)
text(label.pos$x[1]-.02, label.pos$y[1]-.04, labels = "Fish Size (mm)", adj=c(0, NA), srt=302, cex=1.75)

labels <- c("6", "6.5", "7", "7.5", "8")
label.pos <- trans3d((max.x + 0.1), y.axis, min.z, pmat)
text(label.pos$x, label.pos$y, labels=labels, adj=c(0, NA), cex=1.5)
text(label.pos$x[1]+.35, label.pos$y[1]+.025, labels = "Prey Size (mm)", cex=1.75, srt = 16.25)

labels <- c("", as.character(z.axis)[-1])
label.pos <- trans3d(min.x, (min.y - 0.1), z.axis, pmat)
text(label.pos$x, label.pos$y, labels=labels, adj=c(1, NA), cex=1.5)
text(label.pos$x[4]-.06, label.pos$y[4]-.05, labels = "Rel. Selection", cex = 1.75, srt = 270)

#-----------------------------------------------------------------------------#
# End