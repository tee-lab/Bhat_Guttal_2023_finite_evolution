############################
# Plot the steady state solution of the ecological rate modulation example
# introduced in the main text
# To illustrate the fast mechanism of noise-induced selection
# Written by Shikhara Bhat
# JGU Mainz, Germany
# Date: Thu Jan 25 16:08:24 2024
###########################
#Import required packages
library(ggplot2)
library(latex2exp) #to include TeX in ggplot legends

###########################
#Define the functions we will be plotting

pvec <- seq(0.1,0.9,0.05)
#pvec will be used to determine the normalization constant
# and thus the y values of our density plots
#we need it because the actual functions blow up at p=0 and p=1 and are thus
#not normalizable (because p=0 and p=1 are absorbing states and the system
#is consequently not ergodic)

#this is the density function without normalizing
no_freq_dep_nonnormal <- function(p,a,c){
  alpha <- 2*a/c
  return ((1/(p*(1-p)))*exp(-alpha*p))
}

#normalize
no_freq_dep <- function(p,a,c){
  N = sum(sapply(pvec,no_freq_dep_nonnormal,a,c))
  return (no_freq_dep_nonnormal(p,a,c)/N)
}

#this is the function f(x) = 1/(x(1-x))
neutral_exp_nonnormal <- function(p){
  return (1/(p*(1-p)))
}

neutral_exp <- function(p){
  N = sum(sapply(pvec,neutral_exp_nonnormal))
  return (neutral_exp_nonnormal(p)/N)
}

#########################
#plotting

#plot the case where fast mechanism does not operate
p <- ggplot(x=pvec) + geom_function(fun=no_freq_dep,args=list(a=1,c=10),linewidth=1,alpha=1,aes(color='b'))

#plot the case where fast mechanism operates
p <- p + geom_function(fun=no_freq_dep,args=list(a=1,c=0.5),linewidth=1,alpha=1,aes(color='c'))

#plot y = 1/(x(1-x))
p <- p + geom_function(fun=neutral_exp,linewidth=1.2,linetype='dashed',aes(color="a"))

p <- p + scale_color_manual(values=c("#000000","#00aaff","#ff0000"),labels=unname(TeX(c("\\textbf{True neutral expectations (1/p(1-p))}","$|\\kappa / c|$ \\textbf{is small (Fast mechanism is weak/absent)}","$|\\kappa / c|$ \\textbf{is not small (Fast mechanism can operate)}"))))


#aesthetics
p <- p + theme_light() + ylab('Quasi-stationary density') + xlab("Trait frequency (p)")
p <- p + scale_x_continuous(breaks=c(0.05,0.25,0.5,0.75,0.95),limits=c(0.05,0.95))
p <- p + labs(color='')
p <- p + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p + theme(axis.text = element_text(face = 'bold', color = 'black',size = 20))
p <- p + theme(axis.title = element_text(face = 'bold',color = 'black',size = 25))
p <- p + theme(legend.title = element_text(face = 'bold',color = 'black',size = 20))
p <- p + theme(legend.text = element_text(face = 'bold',color = 'black',size = 15))
p <- p + theme(legend.position = c(0.5, 0.8))
p <- p + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',linewidth=0.5))
p <- p + theme(aspect.ratio = 0.8)
p

#save to file
ggsave('./plots/fast_mechanism.svg',p,width=9.42,height=7.26,units='in')