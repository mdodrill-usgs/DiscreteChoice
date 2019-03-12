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
library(fishR)

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
setwd("U:/Desktop/Fish_Git/DiscreteChoice/")
source(paste0(getwd(),"/Code/Functions.R"))

model.set.up(model.name = "Width")
emp.dat = data.frame(taxa = name.key$name2,
                     emp = emp_a_2)
#-----------------------------------------------------------------------------#
# test.1 = bayes_summary(fit.1, par.name = "gamma")

#-----------------------------------------------------------------------------#
name = "gamma"

f1 <- coda::mcmc.list(lapply(1:ncol(fit.1), function(x) coda::mcmc(as.array(fit.1)[,x,])))
f1 <- ggs(f1, family = "gamma")

f1.sub = f1[which(substr(f1$Parameter,1,nchar(name)) == name),]

all = f1.sub

all$id = as.numeric(gsub("[^0-9]", "", as.character(all$Parameter)))
all$taxa = name.key[match(all$id, name.key$num),3]

all$emp = emp.dat[match(all$taxa, emp.dat$taxa),2]

all$ratio = all$value / all$emp

all2 = group_by(all, Parameter) %>%
  summarize(my.mean = mean(ratio),
            upper = quantile(ratio, .975),
            lower = quantile(ratio, .025))


all2$id = as.numeric(gsub("[^0-9]", "", as.character(all2$Parameter)))

all2$taxa = name.key[match(all2$id, name.key$num),2]

#-----------------------------------------------------------------------------#
# Selection Ratios
windows(width = 6.5*1.5, height = 5*1.5, record = T, xpos = 25)

all2$taxa2 = name.key[match(all2$taxa, name.key$name),3]

sub.1 = all2[all2$taxa == "Black Fly Adult",]

sub.2 = all2[all2$taxa != "Black Fly Adult",]


p1 = ggplot(sub.1, aes(y = my.mean, x = taxa2)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  geom_point(position = position_dodge(width = 1), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 1),
                width = 0, size = 1) +
  
  # scale_color_manual(values = c("black", "gray50", "gray")) +
  scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80)) +
  labs(y = "Selection \n Ratio (+ - 95% CRI)", x = "", title = "A")

g1 = p1 + theme(axis.title.x = element_text(size = 18, vjust = -.1),
                axis.title.y = element_text(size = 18, vjust = 4),
                axis.text.x = element_text(size = 16, colour = "black"),
                axis.text.y = element_text(size = 16, colour = "black"),
                title = element_text(size = 18),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect(fill = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                # panel.border = element_rect(colour = "black", fill = NA),
                panel.spacing = unit(.8, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(size = 14, vjust = 1),
                legend.text = element_text(size = 16),
                legend.key = element_rect(fill = "white"),
                legend.title = element_text(size = 16),
                legend.position = "none",
                plot.margin = unit(c(1,0,.5,.75), "lines"))  #T, R, B, L)
g1


p2 = ggplot(sub.2, aes(y = my.mean, x = taxa2)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  geom_point(position = position_dodge(width = 1), size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 1),
                width = 0, size = 1) +
  
  # scale_color_manual(values = c("black", "gray50", "gray")) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8)) +
  labs(y = "", x = "",color = "Fish Size")

g2 = p2 + theme(axis.title.x = element_text(size = 18, vjust = -.1),
                axis.title.y = element_text(size = 18, vjust = 1),
                axis.text.x = element_text(size = 16, colour = "black"),
                axis.text.y = element_text(size = 16, colour = "black"),
                axis.line = element_line(colour = "black"),
                panel.background = element_rect(fill = "white"),
                panel.grid.minor = element_line(colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                # panel.border = element_rect(colour = "black", fill = NA),
                panel.spacing = unit(.8, "lines"),
                strip.background = element_blank(),
                strip.text = element_text(size = 14, vjust = 1),
                legend.text = element_text(size = 16),
                legend.key = element_rect(fill = "white"),
                legend.title = element_text(size = 16),
                legend.position = c(.85, .9),
                plot.margin = unit(c(1,1,.5,0), "lines"))  #T, R, B, L
g2



grid.arrange(g1, g2, ncol = 2, widths = c(1.75/7,6/7))

#-----------------------------------------------------------------------------#
all = f1.sub

all2 = group_by(all, Parameter) %>%
  summarize(my.mean = mean(value),
            upper = quantile(value, .975),
            lower = quantile(value, .025))

all2$id = as.numeric(gsub("[^0-9]", "", as.character(all2$Parameter)))

all2$taxa = name.key[match(all2$id, name.key$num),2]

all2$taxa2 = name.key[match(all2$taxa, name.key$name),3]

emp.dat.3 = emp.dat[order(emp.dat$taxa),]

emp.dat.3$x.st = seq(1:7) - .15
emp.dat.3$x.ed = seq(1:7) + .15

gam_all = all2

#-----------------------------------------------------------------------------#
# diet proportions

p = ggplot(gam_all, aes(y = my.mean * 100, x = taxa2)) +
  geom_point(position = position_dodge(width = 1), size = 3) +
  geom_errorbar(aes(ymin = lower*100, ymax = upper*100),
                position = position_dodge(width = 1),
                width = 0, size = 1) +
  # scale_y_continuous(breaks = c(0,.1,.2,.3,.4,.5,.6,.7), limits = c(0,.725)) +
  # scale_y_continuous(breaks = c(0,.1,.2,.3,.4,.5), limits = c(0,.525)) +
  scale_y_continuous(breaks = c(0,10,20,30,40,50), limits = c(0, 52.5)) +
  geom_segment(data = emp.dat.3, aes(y = emp*100, x = x.st, xend = x.ed, yend = emp*100), size = 3, color = "gray") +
  
  # totally hack to get the points over the black bars, but works;)
  geom_point(position = position_dodge(width = 1), size = 3) +
  geom_errorbar(aes(ymin = lower*100, ymax = upper*100),
                position = position_dodge(width = 1),
                width = 0, size = 1) +
  
  # coord_flip() +
  labs(y = "Predicted Diet \n Percentage (+ - 95% CRI)", x = "", title = "B") #+
# annotate("text", x = .65, y = .725, label = "A", size = 7)
# labs(y = "", x = "", color = "Fish Size")


p3 = p + theme(axis.title.x = element_text(size = 18, vjust = -.1),
               axis.title.y = element_text(size = 18, vjust = 4),
               axis.text.x = element_text(size = 16, colour = "black"),
               axis.text.y = element_text(size = 16, colour = "black"),
               title = element_text(size = 18),
               axis.line = element_line(colour = "black"),
               panel.background = element_rect(fill = "white"),
               panel.grid.minor = element_line(colour = "white"),
               panel.grid.major = element_line(colour = "white"),
               # panel.border = element_rect(colour = "black", fill = NA),
               panel.spacing = unit(.8, "lines"),
               strip.background = element_blank(),
               strip.text = element_text(size = 14, vjust = 1),
               legend.text = element_text(size = 16),
               legend.key = element_rect(fill = "white"),
               legend.title = element_text(size = 16),
               legend.position = c(.875, .9),
               plot.margin = unit(c(1,1,1,.75), "cm")) #t,r,b,l
p3

#-----------------------------------------------------------------------------#
windows(width = 6.5*1.5, height = 5*1.5, record = T, xpos = 25)

p1 = grid.arrange(g1, g2, ncol = 2, widths = c(1.75/7,6/7))


grid.arrange(p1,p3, ncol = 2)




#-----------------------------------------------------------------------------#