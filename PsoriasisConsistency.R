# Methods to Assess Evidence Consistency in Dose-Response Model Based Network Meta-Analysis
# Accompanying R code
# Author: Hugo Pedder

# Whether to save plots or not
saveplots <- TRUE

n.iter <- 100000
n.burnin <- 50000
n.chain <- 4
sdprior <- list(sd="dunif(0,5)")


####### LOAD WORKING VERSION OF MBNMAdose #######
library(MBNMAdose)
library(ggplot2)
library(tidyverse)

#### Prepare Psoriasis dataset (from Warren) ####
library(officer)
library(readr)

# Load psoriasis data from MBNMAdose package
data("psoriasis75")

# Exclude agents with a single dose
psoriasis75 <- psoriasis75[!psoriasis75$agent %in% c("Adalimumab", "Guselkumab"),]

# Exclude studies with a single arm
psori <- psoriasis75 %>% group_by(studyID) %>% filter(n()>1)

# Create MBNMA network object
network <- mbnma.network(psori)





########################
#### Run NMAs ####
########################

set.seed(890421)

#### At the treatment-level ####
plot(network)

# Common treatment effects
nma.f <- nma.run(network, n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain,
                 parameters.to.save=c("d", "resdev", "totresdev"))

# Random treatment effects
nma.r <- nma.run(network, method="random", n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain,
                 parameters.to.save=c("d", "sd", "resdev", "totresdev"),
                 priors=sdprior)


#### At the agent-level ####

netagent <- network
netagent$data.ab$treatment <- netagent$data.ab$agent
netagent$treatments <- netagent$agents

plot(netagent)

# Common treatment effects
nma.agent.f <- nma.run(netagent, n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain,
                       parameters.to.save=c("d", "resdev", "totresdev"))

# Random treatment effects
nma.agent.r <- nma.run(netagent, method="random", n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain,
                       parameters.to.save=c("d", "sd", "resdev", "totresdev"),
                       priors=sdprior)


# Save models
save(nma.f, nma.r, nma.agent.f, nma.agent.r, file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nmaresults.RData")



########################
#### Run MBNMAs ####
########################

set.seed(890421)

# Linear
lin.f <- mbnma.run(network, fun=dpoly(degree=1), n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl",
                   parameters.to.save = c("beta.1", "resdev", "totresdev"))
lin.r <- mbnma.run(network, fun=dpoly(degree=1), method="random", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl",
                   parameters.to.save = c("beta.1", "sd", "resdev", "totresdev"),
                   priors=sdprior)

# Log-linear (not reported in manuscript)
loglin.f <- mbnma.run(network, fun=dloglin(), n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl",
                      parameters.to.save = c("rate", "resdev", "totresdev"))
loglin.r <- mbnma.run(network, fun=dloglin(), method="random", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl",
                      parameters.to.save = c("rate", "sd", "resdev", "totresdev"),
                      priors=sdprior)

# Exponential
exp.f <- mbnma.run(network, fun=dexp(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                   parameters.to.save = c("rate", "resdev", "totresdev"))
exp.r <- mbnma.run(network, fun=dexp(), method="random", pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                   parameters.to.save = c("rate", "sd", "resdev", "totresdev"),
                   priors=sdprior)

# Emax
emax.f <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                    parameters.to.save = c("emax", "ed50", "resdev", "totresdev"))
emax.r <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, method="random",
                    parameters.to.save = c("emax", "ed50", "sd", "resdev", "totresdev"),
                    priors=sdprior)


# Save models
save(lin.f, lin.r, loglin.f, loglin.r, exp.f, exp.r, emax.f, emax.r,
     file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/mbnmaresults.RData")



########################
#### Run NMAs (UME) ####
########################

set.seed(890421)

#### At the treatment-level ####
plot(network)

# Common treatment effects
nma.ume.f <- nma.run(network, n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain, UME = TRUE,
                     parameters.to.save = c("d", "resdev", "totresdev"))

# Random treatment effects
nma.ume.r <- nma.run(network, method="random", n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain, UME=TRUE,
                     parameters.to.save = c("d", "sd", "resdev", "totresdev"),
                     priors=sdprior)


#### At the agent-level ####

netagent <- network
netagent$data.ab$treatment <- netagent$data.ab$agent
netagent$treatments <- netagent$agents

plot(netagent)

# Common treatment effects
nma.agent.ume.f <- nma.run(netagent, n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain, UME=TRUE,
                           parameters.to.save = c("d", "resdev", "totresdev"))

# Random treatment effects
nma.agent.ume.r <- nma.run(netagent, method="random", n.iter=n.iter, n.burnin=n.burnin, pd="pd.kl", n.chain=n.chain, UME=TRUE,
                           parameters.to.save = c("d", "sd", "resdev", "totresdev"),
                           priors=sdprior)


# Save models
save(nma.ume.f, nma.ume.r, nma.agent.ume.f, nma.agent.ume.r,
     file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nmaresults_ume.RData")




########################
#### Run MBNMAs (UME) ####
########################

set.seed(890421)

# Linear
lin.ume.f <- mbnma.run(network, fun=dpoly(degree=1), n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl", UME=TRUE,
                       parameters.to.save = c("beta.1", "resdev", "totresdev"))
lin.ume.r <- mbnma.run(network, fun=dpoly(degree=1), method="random", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl", UME=TRUE,
                       parameters.to.save = c("beta.1", "sd", "resdev", "totresdev"),
                       priors=sdprior)

# Log-linear (not reported in manuscript)
loglin.ume.f <- mbnma.run(network, fun=dloglin(), n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl", UME=TRUE,
                          parameters.to.save = c("rate", "resdev", "totresdev"))
loglin.ume.r <- mbnma.run(network, fun=dloglin(), method="random", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, pd="pd.kl", UME=TRUE,
                          parameters.to.save = c("rate", "sd", "resdev", "totresdev"),
                          priors=sdprior)

# Exponential
exp.ume.f <- mbnma.run(network, fun=dexp(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, UME=TRUE,
                       parameters.to.save = c("rate", "resdev", "totresdev"))
exp.ume.r <- mbnma.run(network, fun=dexp(), method="random", pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, UME=TRUE,
                       parameters.to.save = c("rate", "sd", "resdev", "totresdev"),
                       priors=sdprior)

# Emax
emax.ume.f <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, UME=TRUE,
                        parameters.to.save = c("emax", "ed50", "resdev", "totresdev"))
emax.ume.r <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain, method="random", UME=TRUE,
                        parameters.to.save = c("emax", "ed50", "sd", "resdev", "totresdev"),
                        priors=sdprior)



# Save models
save(lin.ume.f, lin.ume.r, loglin.ume.f, loglin.ume.r, exp.ume.f, exp.ume.r, emax.ume.f, emax.ume.r,
     file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/mbnmaresults_ume.RData")






###############################################################
################# Figures and Tables ################
###############################################################


#### Figure 1: Network plots ####

library(igraph)

# Treatment level
g <- plot(network, v.color = "agent", v.scale = 2, label.distance = 5)

V(g)$name <- sapply(names(V(g)), FUN=function(x) {
  letter <- substring(x,1,2)
  dose <- strsplit(x, split="_", fixed=TRUE)[[1]][2]
  return(paste(letter, dose, sep=" "))
})
V(g)$label.cex <- 0.7
#V(g)$size <- v.size
E(g)$color <- "black"

# Agent level
gagent <- plot(network, v.color = "agent", level="agent", remove.loops = TRUE, v.scale = 2, label.distance = 8)
V(gagent)$label.cex <- 0.7
E(gagent)$color <- "black"
V(gagent)$color <- unique(V(g)$color)
V(gagent)$label.dist <- c(8,8,6,6,8,5,6,9)

legnam <- paste(network$agents, c("(Pl)", "(Br)", "(Et)", "(In)", "(Ix)", "(Se)", "(Ti)", "(Us)"))


if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/Figure1.pdf")
  par(mfrow=c(1,2), mar=c(0, 2.5, 0, 2.5), xpd=TRUE)
  plot(gagent, layout=layout_in_circle(gagent))
  plot(g, layout=layout_in_circle(g))
  legend(x=-1.5, y=-1.7, legend=legnam, pt.bg=unique(V(g)$color),
         pch=21, pt.cex=1, cex=0.8)
  dev.off()
}



cols <- c(RColorBrewer::brewer.pal(6, "Blues")[5], RColorBrewer::brewer.pal(6, "Greens")[5], RColorBrewer::brewer.pal(6, "Reds")[5])
E(g)$lty <- 1
E(g, P=c("Ix 20", "Ix 40"))$color <- cols[1]
E(g, P=c("Ix 20", "Ix 40"))$lty <- 2
E(g, P=c("Et 100", "Us 3.75"))$color <- cols[1]
E(g, P=c("Et 100", "Us 3.75"))$lty <- 2
E(g, P=c("Se 75", "Us 5.62"))$color <- cols[1]
E(g, P=c("Se 75", "Us 5.62"))$lty <- 2
E(g, P=c("Et 100", "Us 7.5"))$color <- cols[1]
E(g, P=c("Et 100", "Us 7.5"))$lty <- 2

E(g, P=c("Ix 30", "Us 5.62"))$color <- cols[3]
E(g, P=c("Ix 30", "Us 5.62"))$lty <- 3

E(g, P=c("Pl 0", "Se 37.5"))$color <- cols[3]
E(g, P=c("Pl 0", "Se 37.5"))$lty <- 3
E(g, P=c("Pl 0", "Ti 8.3"))$color <- cols[3]
E(g, P=c("Pl 0", "Ti 8.3"))$lty <- 3
E(g, P=c("Pl 0", "Ti 16.7"))$color <- cols[3]
E(g, P=c("Pl 0", "Ti 16.7"))$lty <- 3
E(g, P=c("Pl 0", "Br 70"))$color <- cols[3]
E(g, P=c("Pl 0", "Br 70"))$lty <- 3
E(g, P=c("Pl 0", "Br 105"))$color <- cols[3]
E(g, P=c("Pl 0", "Br 105"))$lty <- 3
E(g, P=c("Pl 0", "Et 50"))$color <- cols[3]
E(g, P=c("Pl 0", "Et 50"))$lty <- 3
E(g, P=c("Pl 0", "In 34.5"))$color <- cols[3]
E(g, P=c("Pl 0", "In 34.5"))$lty <- 3
E(g, P=c("Pl 0", "In 57.62"))$color <- cols[3]
E(g, P=c("Pl 0", "In 57.62"))$lty <- 3

E(g)$width <- 3
V(g)$size <- 12


if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/FigureS3.pdf")
  par(mar=c(12, 3.5, 4.1, 4.1))
  plot(g, layout=layout_in_circle(g))
  legend(x=0.8, y=-1.5, legend=legnam, pt.bg=unique(V(g)$color),
         pch=21, pt.cex=1, cex=0.8)
  legend(x=-2, y=-1.5, title=as.expression(bquote(bold("Indirect evidence source"))),
         legend=c("Consistency & dose-response", "Dose-response only"),
         lty=c(2,3), lwd=2, col=cols[c(1,3)], cex=0.8)
  dev.off()
}






###### Table S1: Raw data ######

tab.df <- psori %>% group_by(studyID) %>% mutate(studyindex=cur_group_id())
tab.df$agent <- factor(tab.df$agent, levels=c("Placebo", sort(unique(tab.df$agent[tab.df$agent!="Placebo"]))))
tab.df <- tab.df %>% group_by(studyID) %>% arrange(studyID, agent, dose) %>% mutate(armindex=seq(n()))

tab.df <- data.frame("studyID"=tab.df$studyID,
                     "studyindex"=tab.df$studyindex,
                     "armindex"=tab.df$armindex,
                     "Agent"=tab.df$agent,
                     "Dose"=tab.df$dose,
                     "r"=tab.df$r,
                     "n"=tab.df$N
                     )
write.csv(tab.df, row.names = FALSE, file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMAdoseConsistencyAnalysis/TableA1.csv")



##### Figure S1: Raw data plot ######

psori <- network$data.ab

# Generate additional rows for placebo data so that it is included in each relevant agent panel
psori <- psori %>% group_by(studyID) %>% mutate(temp=ifelse(agent==1, agent[2], agent))
psori <- psori %>% group_by(studyID) %>% mutate(temp2=ifelse(agent==1, agent[length(agent)], agent)) %>%
  mutate(match=ifelse(all(temp==temp2), 1, 0))

newrows <- psori %>% group_by(studyID) %>% slice(1)
newrows <- subset(newrows, match==0 & agent==1)
newrows <-  newrows %>% group_by(studyID) %>% mutate(agent=ifelse(agent==1, temp2, agent))

psori <- psori %>% group_by(studyID) %>% mutate(agent=ifelse(agent==1, temp, agent))
psori <- rbind(psori, newrows)

# Relabel agents
psori$agent <- factor(psori$agent, labels=network$agents[-1])

# Plot raw responses by dose (Figure S1)
g <- ggplot(data=psori, aes(x=dose, y=r/N)) +
  geom_point() +
  facet_wrap(~agent, scales = "free_x") +
  xlab("Dose mg/week") +
  ylab(expression("Proportion of patients with " >= "75% PASI improvement")) +
  theme_mbnma()

if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/FigureS1.pdf")
  plot(g)
  dev.off()
}







########### Figure 2 & 3: Dose-Response Examples #############

# For Figure 2

cols <- RColorBrewer::brewer.pal(n=3, name="Set1")

e0 <- 5
emaxa <- 100
ed50a <- 15
emaxb <- 70
ed50b <- 30

x <- seq(0,50,0.5)
y <- e0 + ((emaxa * x) / (ed50a + x))

dose <- x
theta <- y

y <- e0 + ((emaxb * x) / (ed50b + x))

dose <- c(dose, x)
theta <- c(theta, y)


df.dr <- data.frame("dose"=dose, "theta"=theta, "agent"=c(rep(1,length(x)), rep(2, length(x))))
df.dr$Agent <- factor(df.dr$agent, labels=c("Agent A", "Agent B"))

g <- ggplot(df.dr, aes(x=dose, y=theta, color=Agent, linetype=Agent)) +
  geom_line(size=1) +
  ylim(c(0,90)) +
  scale_color_manual(values=cols[1:2])

# Add points
df.point <- data.frame("dose"=c(0,20,30,0), "theta"=c(e0,
                                                    e0 + ((emaxa * 20) / (ed50a + 20)),
                                                    e0 + ((emaxb * 30) / (ed50b + 30)),
                                                    e0
                                                    ),
                       "Agent"=c("Agent A", "Agent A", "Agent A", "Agent A"))

g <- g + geom_point(data=df.point, color="black", size=3) +
  geom_path(data=df.point, size=1, color="black", linetype="solid") +
  scale_linetype_manual(values=c("dotdash", "dashed")) +
  xlab("Dose") + ylab(expression(theta)) +
  guides(linetype=guide_legend(keywidth = 2.5), color=guide_legend(keywidth = 2.5)) +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14)
        )


# Save plot to powerpoint
if (saveplots==TRUE) {
  library(officer)
  gout <- rvg::dml(code = plot(g))

  officer::read_pptx() %>%
    officer::add_slide() %>%
    officer::ph_with(gout, ph_location_fullsize()) %>%
    base::print(target = "~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/Figure2_temp.pptx")

}




# For Figure 3

g <- ggplot(df.dr, aes(x=dose, y=theta, color=Agent)) +
  geom_line(size=1.3) +
  ylim(c(0,90)) +
  scale_color_manual(values=cols[1:2])


df.point <- data.frame("dose"=c(0,20,45,0,
                                0,30,50,0,
                                45,50),
                       "theta"=c(e0,
                                 e0 + ((emaxa * 20) / (ed50a + 20)),
                                 e0 + ((emaxa * 45) / (ed50a + 45)),
                                 e0,
                                 e0,
                                 e0 + ((emaxb * 30) / (ed50b + 30)),
                                 e0 + ((emaxb * 50) / (ed50b + 50)),
                                 e0,
                                 e0 + ((emaxa * 45) / (ed50a + 45)),
                                 e0 + ((emaxb * 50) / (ed50b + 50))
                                 ),
                       "Study"=c(rep("Study 1", 4), rep("Study 2", 4), rep("Study 3",2)))

g <- g + geom_point(data=df.point, color="black", size=3) +
  geom_path(data=df.point, aes(linetype=Study), size=1, color="black") +
  #scale_linetype_manual(values=c("dotdash", "dashed")) +
  xlab("Dose") + ylab(expression(theta)) +
  guides(linetype=guide_legend(keywidth = 3), color=guide_legend(keywidth = 3)) +
  theme_bw() +
  theme(axis.text=element_text(size=14),
        axis.title.y=element_text(size=16),
        axis.title.x=element_text(size=14),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14)
  )


# Save plot to powerpoint
if (saveplots==TRUE) {
  gout <- rvg::dml(code = plot(g))

  officer::read_pptx() %>%
    officer::add_slide() %>%
    officer::ph_with(gout, ph_location_fullsize()) %>%
    base::print(target = "~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/Figure3_temp.pptx")

}



########## Table 1: Model Fit Results #############

library(xlsx)
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nmaresults.RData")
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/mbnmaresults.RData")
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nmaresults_ume.RData")
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/mbnmaresults_ume.RData")

modfit.df <- data.frame(
  model=c(rep("NMA treatment-level",2),
          rep("NMA agent-level", 2),
          rep("MBNMA linear",2),
          rep("MBNMA exponential",2),
          rep("MBNMA Emax", 2),
          rep("UME treatment-level",2)
          ),
  te=c(rep(c("Common", "Random"), 6)),
  totresdev=round(c(
    nma.f$jagsresult$BUGSoutput$mean$totresdev,
    nma.r$jagsresult$BUGSoutput$mean$totresdev,
    nma.agent.f$jagsresult$BUGSoutput$mean$totresdev,
    nma.agent.r$jagsresult$BUGSoutput$mean$totresdev,
    lin.f$BUGSoutput$mean$totresdev,
    lin.r$BUGSoutput$mean$totresdev,
    exp.f$BUGSoutput$mean$totresdev,
    exp.r$BUGSoutput$mean$totresdev,
    emax.f$BUGSoutput$mean$totresdev,
    emax.r$BUGSoutput$mean$totresdev,
    nma.ume.f$jagsresult$BUGSoutput$mean$totresdev,
    nma.ume.r$jagsresult$BUGSoutput$mean$totresdev
  ),1),
  pd=round(c(
    nma.f$jagsresult$BUGSoutput$pD,
    nma.r$jagsresult$BUGSoutput$pD,
    nma.agent.f$jagsresult$BUGSoutput$pD,
    nma.agent.r$jagsresult$BUGSoutput$pD,
    lin.f$BUGSoutput$pD,
    lin.r$BUGSoutput$pD,
    exp.f$BUGSoutput$pD,
    exp.r$BUGSoutput$pD,
    emax.f$BUGSoutput$pD,
    emax.r$BUGSoutput$pD,
    nma.ume.f$jagsresult$BUGSoutput$pD,
    nma.ume.r$jagsresult$BUGSoutput$pD

  ),1),
  dic=round(c(
    nma.f$jagsresult$BUGSoutput$DIC,
    nma.r$jagsresult$BUGSoutput$DIC,
    nma.agent.f$jagsresult$BUGSoutput$DIC,
    nma.agent.r$jagsresult$BUGSoutput$DIC,
    lin.f$BUGSoutput$DIC,
    lin.r$BUGSoutput$DIC,
    exp.f$BUGSoutput$DIC,
    exp.r$BUGSoutput$DIC,
    emax.f$BUGSoutput$DIC,
    emax.r$BUGSoutput$DIC,
    nma.ume.f$jagsresult$BUGSoutput$DIC,
    nma.ume.r$jagsresult$BUGSoutput$DIC
  ),1),
  sd=c("-", neatCrI(nma.r$jagsresult$BUGSoutput$summary["sd",c(3,5,7)], digits=2),
       "-", neatCrI(nma.agent.r$jagsresult$BUGSoutput$summary["sd",c(3,5,7)], digits=2),
       "-", neatCrI(lin.r$BUGSoutput$summary["sd",c(3,5,7)], digits=2),
       "-", neatCrI(exp.r$BUGSoutput$summary["sd",c(3,5,7)], digits=2),
       "-", "-",
       "-", neatCrI(nma.ume.r$jagsresult$BUGSoutput$summary["sd",c(3,5,7)], digits=2)
       ),
  select=c("Yes", "No",
           "No", "Yes",
           "No", "Yes",
           "No", "Yes",
           "Yes", "No",
           "-", "-")
)

write.xlsx(modfit.df, file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/modelfit.xlsx",
           row.names=FALSE)





###########  Figure 5: Predicted Plots  ###########

load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/mbnmaresults.RData")
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nmaresults.RData")

cols <- RColorBrewer::brewer.pal(5, "Set1")

E0 <- 0.05
n.doses <- 30
pred.emax <- plot(predict(emax.f, E0=E0, n.doses = n.doses), overlay.split = TRUE)
pred.emax <- pred.emax$data
pred.emax <- pred.emax[,c("agent", "dose", "2.5%", "50%", "97.5%")]
pred.emax$model <- "Emax MBNMA"

pred.exp <- plot(predict(exp.r, E0=E0, n.doses = n.doses), overlay.split = FALSE)
pred.exp <- pred.exp$data
pred.exp <- pred.exp[,c("agent", "dose", "2.5%", "50%", "97.5%")]
pred.exp$model <- "Exponential MBNMA"

pred.lin <- plot(predict(lin.r, E0=E0, n.doses = n.doses), overlay.split = FALSE)
pred.lin <- pred.lin$data
pred.lin <- pred.lin[,c("agent", "dose", "2.5%", "50%", "97.5%")]
pred.lin$model <- "Linear MBNMA"

pred.agent <- nma.agent.r$jagsresult$BUGSoutput$sims.matrix[,
                                                            grepl("^d\\[", colnames(nma.agent.r$jagsresult$BUGSoutput$sims.matrix))
                                                            ]
pred.agent <- rescale.link(E0, direction="link", link="logit") + pred.agent
pred.agent <- rescale.link(pred.agent, direction="natural", link="logit")
pred.agent <- t(apply(pred.agent, MARGIN=2, FUN=function(x) {
  quantile(x, probs = c(0.025, 0.5, 0.975))
}))[-1,]
rownames(pred.agent) <- network$agents[-1]

temp <- pred.emax
for (i in seq_along(unique(pred.emax$agent))) {
  agent <- as.character(unique(pred.emax$agent)[i])
  rows <- temp[temp$agent==agent,]
  temp[temp$agent==agent, c("2.5%", "50%", "97.5%")] <- pred.agent[rep(which(rownames(pred.agent)==agent), nrow(rows)),]
}
pred.agent <- temp
pred.agent$model <- "Agent-level NMA"

plot.df <- rbind(pred.agent, pred.exp, pred.lin, pred.emax)
plot.df$Model <- factor(plot.df$model)


g <- ggplot(plot.df, aes(x=dose, y=`50%`, ymin=`2.5%`, ymax=`97.5%`, group=Model, linetype=Model)) +
  geom_ribbon(aes(fill=Model), alpha=0.2) +
  geom_line() +
  facet_wrap(~agent, scales="free_x") +
  scale_fill_manual(values=cols[1:4])



split.df <- nma.f[["jagsresult"]]$BUGSoutput$summary
split.df <- as.data.frame(split.df[grepl("^d\\[[0-9]+\\]",
                                         rownames(split.df)), c(3, 5, 7)])
split.df$treatment <- nma.f[["trt.labs"]]
split.df$agent <- sapply(nma.f[["trt.labs"]], function(x) strsplit(x,
                                                                   split = "_", fixed = TRUE)[[1]][1])
split.df$dose <- as.numeric(sapply(nma.f[["trt.labs"]],
                                   function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2]))
split.df$agent <- factor(split.df$agent, levels = network$agents)
split.df <- split.df[split.df$agent %in% g$data$agent, ]
lE0 <- rescale.link(E0, direction = "link", link = "logit")
split.df[, c(1:3)] <- split.df[, c(1:3)] + lE0
for (i in 1:3) {
  split.df[, i] <- rescale.link(split.df[, i], direction = "natural",
                                link = "logit")
}
split.df$dose[split.df$treatment=="Brodalumab_105"] <- 100
split.df$Model <- NA


g2 <- ggplot(plot.df, aes(x=dose, y=`50%`, ymin=`2.5%`, ymax=`97.5%`, group=Model, linetype=Model)) +
  geom_ribbon(aes(fill=Model), alpha=0.2) +
  geom_line() +
  facet_wrap(~agent, scales="free_x") +
  geom_point(data = split.df, ggplot2::aes(x = dose, y = `50%`), size=1.5) +
  geom_linerange(data = split.df, aes(x = dose, ymin = `2.5%`, ymax = `97.5%`), linetype="solid", size=0.7) +
  xlab("Dose mg/week") + ylab("Proportion of patients with >=75% PASI improvement") +
  theme_mbnma() +
  scale_fill_manual(values=cols[1:4])


if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/Figure5.pdf", width=10)
  plot(g2)
  dev.off()

}








##################################
######### Dev-Dev plots #########
##################################

library(ggplot2)

load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nmaresults.RData")
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nmaresults_ume.RData")
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/mbnmaresults.RData")
load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/mbnmaresults_ume.RData")

# NMA treatment-level

dev.nma <- vector()
dev.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:nma.f$jagsresult$model$data()$NS) {
  narm <- nma.f$jagsresult$model$data()$narm[i]
  dev.nma <- append(dev.nma, nma.f$jagsresult$BUGSoutput$mean$resdev[i,1:narm])
  dev.ume <- append(dev.ume, nma.ume.f$jagsresult$BUGSoutput$mean$resdev[i,1:narm])
  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

trt.df <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.ume, "study"=study, "arm"=arm, "model"="NMA treatment-level")


dev.nma <- vector()
dev.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:nma.f$jagsresult$model$data()$NS) {
  narm <- nma.r$jagsresult$model$data()$narm[i]
  dev.nma <- append(dev.nma, nma.r$jagsresult$BUGSoutput$mean$resdev[i,1:narm])
  dev.ume <- append(dev.ume, nma.ume.r$jagsresult$BUGSoutput$mean$resdev[i,1:narm])
  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

trt.df.r <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.ume, "study"=study, "arm"=arm, "model"="NMA treatment-level (random effects)")

# g <- ggplot(trt.df, aes(x=nma, y  =ume)) +
#    geom_point()



# NMA agent-level (random effects)
netagent <- network
netagent$data.ab$treatment <- netagent$data.ab$agent
netagent$treatments <- netagent$agents

dev.nma <- vector()
dev.ume <- vector()
dev.nma.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:nma.agent.r$jagsresult$model$data()$NS) {
  narm <- nma.agent.r$jagsresult$model$data()$narm[i]
  dev.nma <- append(dev.nma, nma.agent.r$jagsresult$BUGSoutput$mean$resdev[i,1:nma.agent.r$jagsresult$model$data()$narm[i]])
  dev.ume <- append(dev.ume, nma.agent.ume.r$jagsresult$BUGSoutput$mean$resdev[i,1:nma.agent.ume.r$jagsresult$model$data()$narm[i]])
  dev.nma.ume <- append(dev.nma.ume, nma.ume.r$jagsresult$BUGSoutput$mean$resdev[i,1:nma.ume.r$jagsresult$model$data()$narm[i]])
  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

agent.df.r <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm, "model"= "NMA agent-level (random)")

# g <- ggplot(agent.df.r, aes(x=nma, y=ume)) +
#   geom_point()




# NMA agent-level (common effects)
netagent <- network
netagent$data.ab$treatment <- netagent$data.ab$agent
netagent$treatments <- netagent$agents

dev.nma <- vector()
dev.ume <- vector()
dev.nma.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:nma.agent.f$jagsresult$model$data()$NS) {
  narm <- nma.agent.f$jagsresult$model$data()$narm[i]
  dev.nma <- append(dev.nma, nma.agent.f$jagsresult$BUGSoutput$mean$resdev[i,1:nma.agent.f$jagsresult$model$data()$narm[i]])
  dev.ume <- append(dev.ume, nma.agent.ume.f$jagsresult$BUGSoutput$mean$resdev[i,1:nma.agent.ume.f$jagsresult$model$data()$narm[i]])
  dev.nma.ume <- append(dev.nma.ume, nma.ume.f$jagsresult$BUGSoutput$mean$resdev[i,1:nma.ume.f$jagsresult$model$data()$narm[i]])
  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

agent.df.f <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm, "model"= "NMA agent-level (common)")

# g <- ggplot(agent.df.f, aes(x=nma, y=ume)) +
#   geom_point()

agent.df <- rbind(agent.df.r, agent.df.f)





# Linear MBNMA (random effects)

dev.nma <- vector()
dev.ume <- vector()
dev.nma.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:lin.r$model$data()$NS) {
  narm <- lin.r$model$data()$narm[i]

  dev.nma <- append(dev.nma, lin.r$BUGSoutput$mean$resdev[i,1:lin.r$model$data()$narm[i]])
  dev.ume <- append(dev.ume, lin.ume.r$BUGSoutput$mean$resdev[i,1:lin.ume.r$model$data()$narm[i]])
  dev.nma.ume <- append(dev.nma.ume, nma.ume.r$jagsresult$BUGSoutput$mean$resdev[i,1:nma.ume.r$jagsresult$model$data()$narm[i]])

  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

lin.df.r <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm, "model"= "MBNMA Linear (random)")

g <- ggplot(lin.df.r, aes(x=nma, y=ume)) +
  geom_point()



# Linear MBNMA (common effects)

dev.nma <- vector()
dev.ume <- vector()
dev.nma.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:lin.f$model$data()$NS) {
  narm <- lin.f$model$data()$narm[i]

  dev.nma <- append(dev.nma, lin.f$BUGSoutput$mean$resdev[i,1:lin.f$model$data()$narm[i]])
  dev.ume <- append(dev.ume, lin.ume.f$BUGSoutput$mean$resdev[i,1:lin.ume.f$model$data()$narm[i]])
  dev.nma.ume <- append(dev.nma.ume, nma.ume.f$jagsresult$BUGSoutput$mean$resdev[i,1:nma.ume.f$jagsresult$model$data()$narm[i]])

  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

lin.df.f <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm, "model"= "MBNMA Linear (common)")

# g <- ggplot(lin.df.f, aes(x=nma, y=ume)) +
#   geom_point()

lin.df <- rbind(lin.df.r, lin.df.f)




# Exponential MBNMA (random effects)

dev.nma <- vector()
dev.ume <- vector()
dev.nma.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:exp.r$model$data()$NS) {
  narm <- exp.r$model$data()$narm[i]
  dev.nma <- append(dev.nma, exp.r$BUGSoutput$mean$resdev[i,1:exp.r$model$data()$narm[i]])
  dev.ume <- append(dev.ume, exp.ume.r$BUGSoutput$mean$resdev[i,1:exp.ume.r$model$data()$narm[i]])

  dev.nma.ume <- append(dev.nma.ume, nma.ume.r$jagsresult$BUGSoutput$mean$resdev[i,1:nma.ume.r$jagsresult$model$data()$narm[i]])

  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

exp.df <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm, "model"= "MBNMA Exponential (random)")

# g <- ggplot(exp.df, aes(x=nma, y=ume)) +
#   geom_point()



# Exponential MBNMA (common effects)

dev.nma <- vector()
dev.ume <- vector()
dev.nma.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:exp.f$model$data()$NS) {
  narm <- exp.f$model$data()$narm[i]

  dev.nma <- append(dev.nma, exp.f$BUGSoutput$mean$resdev[i,1:exp.f$model$data()$narm[i]])
  dev.ume <- append(dev.ume, exp.ume.f$BUGSoutput$mean$resdev[i,1:exp.ume.f$model$data()$narm[i]])
  dev.nma.ume <- append(dev.nma.ume, nma.ume.f$jagsresult$BUGSoutput$mean$resdev[i,1:nma.ume.f$jagsresult$model$data()$narm[i]])

  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

exp.df.f <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm, "model"= "MBNMA Exponential (common)")

# g <- ggplot(exp.df.f, aes(x=nma, y=ume)) +
#   geom_point()

exp.df <- rbind(exp.df, exp.df.f)




# Emax MBNMA (common effects)

dev.nma <- vector()
dev.ume <- vector()
dev.nma.ume <- vector()
study <- vector()
arm <- vector()
for (i in 1:emax.f$model$data()$NS) {
  narm <- emax.f$model$data()$narm[i]
  dev.nma <- append(dev.nma, emax.f$BUGSoutput$mean$resdev[i,1:emax.f$model$data()$narm[i]])
  dev.ume <- append(dev.ume, emax.ume.f$BUGSoutput$mean$resdev[i,1:emax.ume.f$model$data()$narm[i]])
  dev.nma.ume <- append(dev.nma.ume, nma.ume.f$jagsresult$BUGSoutput$mean$resdev[i,1:nma.ume.f$jagsresult$model$data()$narm[i]])
  study <- append(study, rep(i, narm))
  arm <- append(arm, 1:narm)
}

emax.df <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm, "model"= "MBNMA Emax")

# g <- ggplot(emax.df, aes(x=nma, y=ume)) +
#   geom_point()



##### Plot all ####

dev.df <- rbind(trt.df, trt.df.r, agent.df, lin.df, exp.df, emax.df)

modlevs <- c("MBNMA Emax", "MBNMA Exponential (common)", "MBNMA Exponential (random)",
             "MBNMA Linear (common)", "MBNMA Linear (random)", "NMA agent-level (common)", "NMA agent-level (random)",
             "NMA treatment-level", "NMA treatment-level (random effects)"
             )

modlabs <- c("MBNMA Emax (common effects)",
          "MBNMA Exponential (common effects)",
          "MBNMA Exponential (random effects)",
          "MBNMA Linear (common effects)",
          "MBNMA Linear (random effects)",
          "NMA agent-level (common effects)",
          "NMA agent-level (random effects)",
          "NMA treatment-level (common effects)",
          "NMA treatment-level (random effects)"
)
dev.df$model <- factor(dev.df$model, labels=modlabs, levels=modlevs)


dev.df$split <- dev.df$nma[dev.df$model=="NMA treatment-level (common effects)"]
dev.df$split[grepl("random", dev.df$model)] <- rep(dev.df$nma[dev.df$model=="NMA treatment-level (random effects)"], 4)

dev.df$change <- dev.df$split - dev.df$nma.ume

gmbnma <- ggplot(dev.df, aes(y=ume, x=nma)) +
  geom_point() +
  facet_wrap(~ factor(model), scales = "fixed") +
  xlab("Consistency model") + ylab("MBNMA/NMA UME model") +
  xlim(c(0,10)) + ylim(c(0,10))

gume <- ggplot(dev.df, aes(y=nma.ume, x=nma)) +
  geom_point() +
  facet_wrap(~ factor(model), scales = "fixed") +
  xlab("Consistency model") + ylab("NMA UME model") +
  xlim(c(0,10)) + ylim(c(0,10))

gnma <- ggplot(dev.df, aes(y=split, x=nma)) +
  geom_point() +
  facet_wrap(~ factor(model), scales = "fixed") +
  xlab("Consistency model") + ylab("NMA treatment-level model") +
  xlim(c(0,10)) + ylim(c(0,10))




########### Figure 6: Dev-dev plots ##########

dev.df$incon <- 0
dev.df$incon[dev.df$study==16 & dev.df$arm==1] <- 1
dev.df$incon[dev.df$study==15 & dev.df$arm==1] <- 1
dev.df$incon <- factor(dev.df$incon)

mods <- c("MBNMA Emax (common effects)",
          "MBNMA Exponential (common effects)",
          "NMA agent-level (common effects)",
          "NMA treatment-level (common effects)")

g1 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=split, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: NMA treatment-level model") +
  xlim(c(0,15)) + ylim(c(0,15)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_mbnma() +
  theme(legend.position = "none")


dev.df$incon <- 0
dev.df$incon[dev.df$study==16 & dev.df$arm==1] <- 1
dev.df$incon[dev.df$study==15 & dev.df$arm==1] <- 1
dev.df$incon <- factor(dev.df$incon)

g2 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=nma.ume, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: UME treatment-level model") +
  xlim(c(0,15)) + ylim(c(0,15)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_mbnma() +
  theme(legend.position = "none")


if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/Figure6.pdf", height=10)
  gridExtra::grid.arrange(g1, g2, nrow=2)
  dev.off()
}



###### Figure S2: Dev-dev plots ######

dev.df$incon <- 0
dev.df$incon[dev.df$study==16 & dev.df$arm==1] <- 1
dev.df$incon[dev.df$study==15 & dev.df$arm==1] <- 1
dev.df$incon <- factor(dev.df$incon)

mods <- c("MBNMA Emax (common effects)",
          "MBNMA Exponential (random effects)",
          "NMA agent-level (random effects)",
          "NMA treatment-level (common effects)")

g1 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=split, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: NMA treatment-level model") +
  xlim(c(0,5)) + ylim(c(0,5)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_mbnma() +
  theme(legend.position = "none")

dev.df$incon <- 0
dev.df$incon[dev.df$study==16 & dev.df$arm==1] <- 1
dev.df$incon[dev.df$study==15 & dev.df$arm==1] <- 1
dev.df$incon <- factor(dev.df$incon)

g2 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=nma.ume, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: UME treatment-level model") +
  xlim(c(0,5)) + ylim(c(0,5)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_mbnma() +
  theme(legend.position = "none")


if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/FigureS2.pdf", height=10)
  gridExtra::grid.arrange(g1, g2, nrow=2)
  cowplot::plot_grid(g1, g2, nrow=2, labels=c("a)", "b)"), label_size = 20,
                     label_fontfamily = "serif",
                     label_fontface = "plain")
  dev.off()

}






########## Figure S3: Linear dev-dev plots #########

mods <- c("MBNMA Linear (common effects)",
          "MBNMA Linear (random effects)"
          )

g1 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=split, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: NMA treatment-level model") +
  xlim(c(0,5)) + ylim(c(0,5)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_mbnma() +
  theme(legend.position = "none")


dev.df$incon <- 0
dev.df$incon[dev.df$study==16 & dev.df$arm==1] <- 1
dev.df$incon[dev.df$study==15 & dev.df$arm==1] <- 1
dev.df$incon <- factor(dev.df$incon)

g2 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=nma.ume, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: UME treatment-level model") +
  xlim(c(0,5)) + ylim(c(0,5)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_mbnma() +
  theme(legend.position = "none")


if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/FigureSX.pdf", height=7)
  gridExtra::grid.arrange(g1, g2, nrow=2)
  dev.off()
}






#######################################
######## NODE-SPLITTING ########
#######################################

library(forestplot)
library(dplyr)

######### Manuscript-specific ode-split functions ##########

checkdr <- function(comps, nparams=1) {
  paths <- comps$path
  index <- grep("drparams", paths)
  index <- index[sapply(paths[index], function(x) ifelse(all(strsplit(x, split=" ")[[1]][2:3] >=nparams), TRUE, FALSE))]

  indices <- c(which(!grepl("drparams", paths)), index)

  return(comps[sort(indices),])
}

forestfun <- function(cons.ns, labs=NULL, cex=1) {

  nanal <- length(unique(cons.ns$model))

  cons.ns$Evidence[cons.ns$Evidence %in% c("MBNMA", "NMA")] <- "Overall"

  cons.ns <- arrange(cons.ns, Comparison, model, Evidence)

  names(cons.ns)[4:5] <- c("l95", "u95")

  if (!is.null(labs)) {
    cons.ns$Comparison <- as.character(factor(cons.ns$Comparison, labels=labs))
  }

  cons.ns$p.value <- format(cons.ns$p.value)
  cons.ns$p.value[cons.ns$p.value=="0.000"] <- "<0.001"


  # Remove surplus labels
  row <- nrow(cons.ns)
  cons.ns$model[c(1:row)[!c(1:row) %in% seq(2,row,3)]] <- NA
  cons.ns$p.value[c(1:row)[!c(1:row) %in% seq(2,row,3)]] <- NA
  cons.ns$Comparison[c(1:row)[!c(1:row) %in% seq(1,row, 3*nanal)]] <- NA


  # Add blank between groups
  temp <- cons.ns[0,]
  for (i in 1:(nrow(cons.ns)/3)) {
    seg <- rbind(cons.ns[((i-1)*3+1):(i*3),], rep(NA,7))
    temp <- rbind(temp, seg)
  }
  cons.ns <- temp

  # Add extra blank between comps
  comprow <- which(!is.na(cons.ns$Comparison))
  if (length(comprow)>1) {
    temp <- cons.ns[0,]
    for (i in 1:(length(comprow)-1)) {
      temp <- rbind(temp, rep(NA,7), cons.ns[comprow[i]:(comprow[i+1]-1),])
    }
    temp <- rbind(temp, rep(NA,7), cons.ns[comprow[length(comprow)]:nrow(cons.ns),])
    cons.ns <- temp

    # Summary indicator
    sumind <- rep(F, nrow(cons.ns))
    sumind[which(!is.na(cons.ns$Comparison))-2] <- T
  } else {

    # Summary indicator
    sumind <- rep(F, nrow(cons.ns))
  }



  forestplot(labeltext=cbind(cons.ns$Comparison, cons.ns$model, cons.ns$Evidence, cons.ns$p.value),
             mean=cons.ns$Median, lower=cons.ns$l95, upper=cons.ns$u95,
             boxsize=0.2, graph.pos=4,
             xlab="Log-Odds Ratio (95% CrI)", hrzl_lines = TRUE, is.summary=sumind,
             txt_gp = fpTxtGp(label=gpar(cex=cex), xlab=gpar(cex=cex), ticks=gpar(cex=cex))
  )
}



######## Run node-split analyses ########

comps <- inconsistency.loops(network$data.ab, incldr = TRUE)

# Nodesplit NMA models (treatment-level) with common and random treatment effects
split.nma <- nma.nodesplit(network, n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain)
split.nma.r <- nma.nodesplit(network, method="random", priors=sdprior,
                             n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain)


# Nodesplit linear MBNMA with common and random treatment effects
split.lin <- mbnma.nodesplit(network, fun=dpoly(degree = 1), comparisons = checkdr(comps, 2),
                             n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain)
split.lin.r <- mbnma.nodesplit(network, fun=dpoly(degree = 1), method="random", comparisons = checkdr(comps, 2),
                               n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                               priors=sdprior)

# Nodesplit loglinear MBNMA with common and random treatment effects (not included in manuscript)
split.loglin <- mbnma.nodesplit(network, fun=dloglin(), n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                             comparisons = checkdr(comps, 2))
split.loglin.r <- mbnma.nodesplit(network, fun=dloglin(), method="random", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                               comparisons = checkdr(comps, 2), priors=sdprior)


# Nodesplit exponential MBNMA with common and random treatment effects
split.exp <- mbnma.nodesplit(network, fun=dexp(), n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                             comparisons = checkdr(comps, 2))
split.exp.r <- mbnma.nodesplit(network, fun=dexp(), method="random", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                               comparisons = checkdr(comps, 2), priors=sdprior)


# Nodesplit Emax MBNMA with common treatment effects
split.emax <- mbnma.nodesplit(network, fun=demax(), n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                              comparisons = checkdr(comps, 3))


# Save objects
save(split.nma, split.nma.r,
     split.lin, split.lin.r,
     split.loglin, split.loglin.r,
     split.exp, split.exp.r,
     split.emax,
     file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nodesplits.Rdata")




###### Figures 7 and S5: Node-split forest plots #######

load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nodesplits.Rdata")


############## Figure S7 ################

nma.ns <- summary.nodesplit(split.nma)
nma.ns$model <- "Treatment-level NMA"

exp.ns <- summary.nodesplit(split.exp.r)
exp.ns$model <- "Exponential MBNMA"

emax.ns <- summary.nodesplit(split.emax)
emax.ns$model <- "Emax MBNMA"

cons.ns.1 <- rbind(nma.ns,
                 emax.ns[emax.ns$Comparison %in% c(nma.ns$Comparison),],
                 exp.ns[exp.ns$Comparison %in% c(nma.ns$Comparison),]
)


nma.ns <- summary.nodesplit(split.nma)
nma.ns$model <- "Treatment-level NMA"

exp.ns <- summary.nodesplit(split.exp.r)
exp.ns$model <- "Exponential MBNMA"

emax.ns <- summary.nodesplit(split.emax)
emax.ns$model <- "Emax MBNMA"

cons.ns <- emax.ns[!emax.ns$Comparison %in% c(nma.ns$Comparison),]

cons.ns.2 <- rbind(cons.ns,
                 exp.ns[exp.ns$Comparison %in% c(cons.ns$Comparison),]
)


nma.ns <- summary.nodesplit(split.nma)
nma.ns$model <- "Treatment-level NMA"

exp.ns <- summary.nodesplit(split.exp.r)
exp.ns$model <- "Exponential MBNMA"

emax.ns <- summary.nodesplit(split.emax)
emax.ns$model <- "Emax MBNMA"

cons.ns.3 <- exp.ns[!exp.ns$Comparison %in% c(emax.ns$Comparison),]


cons.ns <- rbind(cons.ns.1, cons.ns.2, cons.ns.3)



# Start plotting procedure

cons.ns$Evidence[cons.ns$Evidence %in% c("MBNMA", "NMA")] <- "Overall"

cons.ns <- cons.ns %>% group_by(Comparison) %>% mutate(count=n())

cons.ns <- arrange(cons.ns, desc(count), Comparison, model, Evidence)

names(cons.ns)[4:5] <- c("l95", "u95")

# Write labels
labs <- unique(cons.ns$Comparison)
labs <- lapply(labs, function(x) {
  trt1 <- strsplit(x, split="_")[[1]][1]
  trt1 <- substr(trt1, 1, 2)

  dose2 <- strsplit(x, split="_")[[1]][3]

  y <- strsplit(x, split="_")[[1]][2]
  dose1 <- strsplit(y, split=" ")[[1]][1]
  trt2 <- strsplit(y, split=" ")[[1]][3]
  trt2 <- substr(trt2, 1, 2)

  return(paste0(trt1, " ", dose1, " vs ", trt2, " ", dose2))
  })

for (i in seq_along(labs)){
  labs[i] <- paste0("#", i, ": ", labs[i])
}

cons.ns$Comparison <- as.character(factor(cons.ns$Comparison, levels=unique(cons.ns$Comparison), labels=c(
  as.expression(labs[1]),
  as.expression(labs[2]),
  as.expression(labs[3]),
  as.expression(labs[4]),

  as.expression(labs[5]),

  as.expression(labs[6]),
  as.expression(labs[7]),
  as.expression(labs[8]),
  as.expression(labs[9]),
  as.expression(labs[10]),
  as.expression(labs[11]),
  as.expression(labs[12]),
  as.expression(labs[13])
)))

cons.ns$p.value <- format(cons.ns$p.value)
cons.ns$p.value[cons.ns$p.value=="0.000"] <- "<0.001"


# Remove surplus labels
cons.ns <- cons.ns %>% group_by(Comparison, model) %>% mutate(count=seq(n()))
cons.ns$model[cons.ns$count!=2] <- NA
cons.ns$p.value[cons.ns$count!=2] <- NA

cons.ns <- cons.ns %>% group_by(Comparison) %>% mutate(count=seq(n()))
cons.ns$Comparison[cons.ns$count>1] <- NA


# Add blank between models
index <- which(cons.ns$Evidence=="Overall")
temp <- cons.ns[1:(index[1]),]
for (i in 1:(length(index)-1)) {
  print(i)
  temp[nrow(temp)+1,] <- NA
  temp <- rbind(temp, cons.ns[(index[i]+1):index[i+1],])
}
cons.ns <- temp


# Add blank between Comparisons
index <- which(!is.na(cons.ns$Comparison))
index <- c(index, nrow(cons.ns)+1)
temp <- cons.ns[0,]
index <- index[]
for (i in 1:(length(index)-1)) {
  print(i)
  temp[nrow(temp)+1,] <- NA
  temp <- rbind(temp, cons.ns[(index[i]):(index[i+1]-1),])
}
temp[nrow(temp)+1,] <- NA
cons.ns <- temp


# Add dividing lines
index <- which(!is.na(cons.ns$Comparison))[-1]
sumind <- rep(FALSE, nrow(cons.ns))
sumind[index-2] <- TRUE


if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/Figure7.pdf", height=11)
  cex <- 0.6
  forestplot(labeltext=cbind(c("Comparison", rep(NA, 1), cons.ns$Comparison),
                             c("Model", rep(NA, 1), cons.ns$model),
                             c("Evidence", rep(NA, 1), cons.ns$Evidence),
                             c("p-value", rep(NA, 1), cons.ns$p.value)),
             mean=c(rep(NA,2),cons.ns$Median), lower=c(rep(NA,2), cons.ns$l95), upper=c(rep(NA,2), cons.ns$u95),
             boxsize=0.4, graph.pos=4,
             xlab="Log-Odds Ratio (95% CrI)", hrzl_lines = TRUE, is.summary=c(TRUE, FALSE, sumind),
             lwd.ci	= 1.7,
             col = fpColors(all.elements = "black"),
             #lineheight=unit(0.1, "cm")
             txt_gp = fpTxtGp(label=gpar(cex=cex), xlab=gpar(cex=cex), ticks=gpar(cex=cex))
  )
  dev.off()
}




############# Figure S5: Supplementary Forest including linear ################

load(file="~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/MBNMADoseConsistencyAnalysis/Model Results/nodesplits.Rdata")

nma.ns <- summary.nodesplit(split.nma)
nma.ns$model <- "Treatment-level NMA"

lin.ns <- summary.nodesplit(split.lin.r)
lin.ns$model <- "Linear MBNMA"

exp.ns <- summary.nodesplit(split.exp.r)
exp.ns$model <- "Exponential MBNMA"

emax.ns <- summary.nodesplit(split.emax)
emax.ns$model <- "Emax MBNMA"

cons.ns.1 <- rbind(nma.ns,
                   lin.ns[lin.ns$Comparison %in% c(nma.ns$Comparison),],
                   emax.ns[emax.ns$Comparison %in% c(nma.ns$Comparison),],
                   exp.ns[exp.ns$Comparison %in% c(nma.ns$Comparison),]
)


cons.ns <- emax.ns[!emax.ns$Comparison %in% c(nma.ns$Comparison),]

cons.ns.2 <- rbind(cons.ns,
                   exp.ns[exp.ns$Comparison %in% c(cons.ns$Comparison),]
)


cons.ns <- exp.ns[!exp.ns$Comparison %in% c(emax.ns$Comparison),]
cons.ns.3 <- rbind(cons.ns,
                   lin.ns[!lin.ns$Comparison %in% c(emax.ns$Comparison),])


cons.ns <- rbind(cons.ns.1, cons.ns.2, cons.ns.3)



# Start plotting procedure

cons.ns$Evidence[cons.ns$Evidence %in% c("MBNMA", "NMA")] <- "Overall"

cons.ns <- cons.ns %>% group_by(Comparison) %>% mutate(count=n())

cons.ns <- arrange(cons.ns, desc(count), Comparison, model, Evidence)

names(cons.ns)[4:5] <- c("l95", "u95")


# Write labels
labs <- unique(cons.ns$Comparison)
labs <- lapply(labs, function(x) {
  trt1 <- strsplit(x, split="_")[[1]][1]
  trt1 <- substr(trt1, 1, 2)

  dose2 <- strsplit(x, split="_")[[1]][3]

  y <- strsplit(x, split="_")[[1]][2]
  dose1 <- strsplit(y, split=" ")[[1]][1]
  trt2 <- strsplit(y, split=" ")[[1]][3]
  trt2 <- substr(trt2, 1, 2)

  return(paste0(trt1, " ", dose1, " vs ", trt2, " ", dose2))
})

for (i in seq_along(labs)){
  labs[i] <- paste0("#", i, ": ", labs[i])
}

cons.ns$Comparison <- as.character(factor(cons.ns$Comparison, levels=unique(cons.ns$Comparison), labels=c(
  as.expression(labs[1]),
  as.expression(labs[2]),
  as.expression(labs[3]),
  as.expression(labs[4]),

  as.expression(labs[5]),

  as.expression(labs[6]),
  as.expression(labs[7]),
  as.expression(labs[8]),
  as.expression(labs[9]),
  as.expression(labs[10]),
  as.expression(labs[11]),
  as.expression(labs[12]),
  as.expression(labs[13])
)))

cons.ns$p.value <- format(cons.ns$p.value)
cons.ns$p.value[cons.ns$p.value=="0.000"] <- "<0.001"


# Remove surplus labels
cons.ns <- cons.ns %>% group_by(Comparison, model) %>% mutate(count=seq(n()))
cons.ns$model[cons.ns$count!=2] <- NA
cons.ns$p.value[cons.ns$count!=2] <- NA

cons.ns <- cons.ns %>% group_by(Comparison) %>% mutate(count=seq(n()))
cons.ns$Comparison[cons.ns$count>1] <- NA


# Add blank between models
index <- which(cons.ns$Evidence=="Overall")
temp <- cons.ns[1:(index[1]),]
for (i in 1:(length(index)-1)) {
  print(i)
  temp[nrow(temp)+1,] <- NA
  temp <- rbind(temp, cons.ns[(index[i]+1):index[i+1],])
}
cons.ns <- temp


# Add blank between Comparisons
index <- which(!is.na(cons.ns$Comparison))
index <- c(index, nrow(cons.ns)+1)
temp <- cons.ns[0,]
index <- index[]
for (i in 1:(length(index)-1)) {
  print(i)
  temp[nrow(temp)+1,] <- NA
  temp <- rbind(temp, cons.ns[(index[i]):(index[i+1]-1),])
}
temp[nrow(temp)+1,] <- NA
cons.ns <- temp


# Add dividing lines
index <- which(!is.na(cons.ns$Comparison))[-1]
sumind <- rep(FALSE, nrow(cons.ns))
sumind[index-2] <- TRUE



if (saveplots==TRUE) {
  pdf("~/MBNMA/MBNMA Manuscripts/Dose-Response Consistency/Graphs/FigureS4.pdf", height=16)
  cex <- 0.6
  forestplot(labeltext=cbind(c("Comparison", rep(NA, 1), cons.ns$Comparison),
                             c("Model", rep(NA, 1), cons.ns$model),
                             c("Evidence", rep(NA, 1), cons.ns$Evidence),
                             c("p-value", rep(NA, 1), cons.ns$p.value)),
             mean=c(rep(NA,2),cons.ns$Median), lower=c(rep(NA,2), cons.ns$l95), upper=c(rep(NA,2), cons.ns$u95),
             boxsize=0.4, graph.pos=4,
             xlab="Log-Odds Ratio (95% CrI)", hrzl_lines = TRUE, is.summary=c(TRUE, FALSE, sumind),
             lwd.ci	= 1.7,
             col = fpColors(all.elements = "black"),
             txt_gp = fpTxtGp(label=gpar(cex=cex), xlab=gpar(cex=cex), ticks=gpar(cex=cex))
  )
  dev.off()
}









#################################################################
############ SENSITIVITY ANALYSES ###############
#################################################################

###### Sensitivity of Emax to Wishart prior ######

emax.f.wish1 <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                          parameters.to.save = c("emax", "ed50"),
                          omega=matrix(c(16,9.6,9.6,9), nrow=2))


emax.f.wish2 <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                    parameters.to.save = c("emax", "ed50"), priors=list(inv.R="dwish(omega[,],20)"))


# Implies correlation of 0.8 between Emax and ED50
emax.f.wish3 <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                          parameters.to.save = c("emax", "ed50"),
                          omega=matrix(c(16,9.6,9.6,9), nrow=2), priors=list(inv.R="dwish(omega[,],20)"))


# Emax model with no correlation between dose-response parameters Does not converge
emax.f.nowish <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                          parameters.to.save = c("emax", "ed50"), cor=FALSE)



######### Figure S6: Sensitivity of results to Wishart prior specification ############

modlist <- list("base"=emax.f, "wishart1"=emax.f.wish1, "wishart2"=emax.f.wish2, "wishart3"=emax.f.wish3)
model <- vector()
agent <- vector()
res <- matrix(ncol=3)
for (i in seq_along(modlist)) {
  ed50 <- modlist[[i]]$BUGSoutput$summary[grepl("^ed50", rownames(modlist[[i]]$BUGSoutput$summary)),]
  emax <- modlist[[i]]$BUGSoutput$summary[grepl("^emax", rownames(modlist[[i]]$BUGSoutput$summary)),]

  agent <- append(agent, rep(modlist[[i]]$network$agents[-1],2))
  model <- append(model, rep(names(modlist)[i], nrow(ed50) + nrow(emax)))

  res <- rbind(res, ed50[,c(3,5,7)])
  res <- rbind(res, emax[,c(3,5,7)])
}
res <- res[-1,]
sens.df <- data.frame(model=model, agent=agent, param=rownames(res))
sens.df <- cbind(sens.df, res)

sens.df$param <- gsub("\\[[0-9]\\]", "", sens.df$param)
rownames(sens.df) <- NULL

sens.df$Priors <- factor(sens.df$model, labels=c("Base-case (no correlation, k=2)",
                                                 "Wishart 1 (strong correlation, k=2)",
                                                 "Wishart 2 (no correlation, k=20)",
                                                 "Wishart 3 (strong correlation, k=20)"))
sens.df$param <- factor(sens.df$param, labels=c("ED50", "Emax"))
sens.df$agent <- factor(sens.df$agent)

dodge <- position_dodge(width=0.5)
cols <- RColorBrewer::brewer.pal(5, "Set1")
ggplot(sens.df, aes(x=`50%`, y=agent, color=Priors, linetype=Priors, shape=Priors)) +
  geom_point(position=dodge) +
  geom_errorbar(aes(xmin=`2.5%`, xmax=`97.5%`), position=dodge, width=0.2) +
  facet_wrap(~param, scales = "free_x") +
  theme_mbnma() +
  scale_color_manual(values=cols[1:n_distinct(sens.df$model)]) +
  xlab("Effect size") + ylab("Agent")



###### Figure S7: Sensitivity of Exponential random effects model to between-study SD prior specification #######

exp.r.sd <- mbnma.run(network, fun=dexp(), method="random", pd="pd.kl", n.iter=n.iter, n.burnin=n.burnin, n.chain=n.chain,
                      parameters.to.save = c("rate", "sd", "resdev", "totresdev"),
                      priors=list(sd="dnorm(0,0.1) T(0,)"))


modlist <- list("base"=exp.r, "Half-normal"=exp.r.sd)
model <- vector()
agent <- vector()
res <- matrix(ncol=3)
for (i in seq_along(modlist)) {
  lambda <- modlist[[i]]$BUGSoutput$summary[grepl("^rate", rownames(modlist[[i]]$BUGSoutput$summary)),]

  agent <- append(agent, modlist[[i]]$network$agents[-1])
  model <- append(model, rep(names(modlist)[i], nrow(lambda)))

  res <- rbind(res, lambda[,c(3,5,7)])
}
res <- res[-1,]
sens.df <- data.frame(model=model, agent=agent, param=rownames(res))
sens.df <- cbind(sens.df, res)

sens.df$param <- gsub("\\[[0-9]\\]", "", sens.df$param)
rownames(sens.df) <- NULL

sens.df$Priors <- factor(sens.df$model, labels=c("Base-case", "Half-Normal"))
sens.df$param <- factor(sens.df$param, labels=c("Lambda"))
sens.df$agent <- factor(sens.df$agent)

dodge <- position_dodge(width=0.5)
cols <- RColorBrewer::brewer.pal(5, "Set1")
ggplot(sens.df, aes(x=`50%`, y=agent, color=Priors, linetype=Priors, shape=Priors)) +
  geom_point(position=dodge) +
  geom_errorbar(aes(xmin=`2.5%`, xmax=`97.5%`), position=dodge, width=0.2) +
  facet_wrap(~param) +
  theme_mbnma() +
  scale_color_manual(values=cols[1:n_distinct(sens.df$model)]) +
  xlab("Effect size") + ylab("Agent")





