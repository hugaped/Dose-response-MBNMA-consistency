# Methods to Assess Evidence Consistency in Dose-Response Model Based Network Meta-Analysis
# Accompnaying R code
# Author: Hugo Pedder

# Requires devtools to be installed
devtools::install_github("hugaped/MBNMAdose", ref="dev_version")
library(MBNMAdose)
library(tidyverse)

# Load psoriasis data
data("psoriasis75")

# Exclude agents with a single dose
psoriasis75 <- psoriasis75[!psoriasis75$agent %in% c("Adalimumab", "Guselkumab"),]

# Exclude studies with a single arm
psoriasis75 <- psoriasis75 %>% group_by(studyID) %>% filter(n()>1)

# Create MBNMA network object
network <- mbnma.network(psoriasis75)




############## Network plots ################

# At treatment level
plot(network, v.color = "agent", label.distance = 5, v.scale = 2)

# At agent level
plot(network, level="agent", v.color = "agent", label.distance=5, v.scale=2, remove.loops = TRUE)



########## Plot raw data by dose (Figure S1) ##########

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
  theme_bw()

plot(g)




###########################################################
################## Run Consistency Models #################
###########################################################

################## Run NMA Models #################

# Common effects "split" NMA at treatment-level
nma.c <- nma.run(network, n.iter=15000, pd="pd.kl", n.chain=4,
                 parameters.to.save=c("d", "totresdev", "resdev"))

# Random effects "split" NMA at treatment-level
nma.r <- nma.run(network, method="random", n.iter=50000, n.burnin=30000, pd="pd.kl", n.chain=4,
                 parameters.to.save=c("d", "totresdev", "resdev", "sd"))


# "Lump" doses of agents together
netagent <- network
netagent$data.ab$treatment <- netagent$data.ab$agent
netagent$treatments <- netagent$agents

# Common effects "lumped" NMA at agent-level
nma.agent.c <- nma.run(netagent, n.iter=50000, n.burnin=30000, pd="pd.kl", n.chain=4,
                       parameters.to.save=c("d", "totresdev", "resdev"))

# Random effects "lumped" NMA at agent-level
nma.agent.r <- nma.run(netagent, method="random", n.iter=50000, n.burnin=30000, pd="pd.kl", n.chain=4,
                       parameters.to.save=c("d", "totresdev", "resdev", "sd"))



############# Run MBNMA Models #############

# Common effects exponential MBNMA
exp.c <- mbnma.run(network, fun=dexp(), pd="pd.kl", n.iter=50000, n.burnin=30000, n.chain=4,
                           parameters.to.save = c("rate", "totresdev", "resdev"))

# Random effects exponential MBNMA
exp.r <- mbnma.run(network, fun=dexp(), method="random", pd="pd.kl", n.iter=50000, n.burnin=30000, n.chain=4,
                           parameters.to.save = c("rate", "totresdev", "resdev", "sd"))

# Common effects Emax MBNMA
emax.c <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=100000, n.burnin=50000, n.chain=4,
                     parameters.to.save = c("emax", "ed50", "totresdev", "resdev"))

# Random effects Emax MBNMA (no convergence)
emax.r <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=100000, n.burnin=50000, n.chain=4, method="random",
                     parameters.to.save = c("emax", "ed50", "totresdev", "resdev", "sd"))




########################################################
############### Plot Predicted Responses ###############
########################################################

E0 <- 0.05 # Define E0 (% responders on placebo / dose=0)
n.doses <- 30 # Set number of doses for which to make predictions

# Simple predictions from Emax MBNMA with split treatment-level NMA predictions overlaid
pred <- predict(emax.c, E0=E0, n.doses = n.doses)
plot(pred, overlay.split=TRUE)



####### Generate predictions for multiple models (Figure 5) #######

# Predictions for Emax MBNMA
pred.emax <- plot(predict(emax.c, E0=E0, n.doses = n.doses), overlay.split = TRUE)
pred.emax <- pred.emax$data
pred.emax <- pred.emax[,c("agent", "dose", "2.5%", "50%", "97.5%")]
pred.emax$model <- "Emax MBNMA"

# Predictions for exponential MBNMA
pred.exp <- plot(predict(exp.r, E0=E0, n.doses = n.doses), overlay.split = FALSE)
pred.exp <- pred.exp$data
pred.exp <- pred.exp[,c("agent", "dose", "2.5%", "50%", "97.5%")]
pred.exp$model <- "Exponential MBNMA"

# Predictions for agent-level "lumped" NMA (must be manually calcualted from NMA model)
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

# Bind predictions from MBNMA and agent-level models together into a single data frame
plot.df <- rbind(pred.agent, pred.exp, pred.emax)
plot.df$Model <- factor(plot.df$model)


# Generate predictions for "split" treatment-level NMA (must be manually calculated from NMA model)
split.df <- nma.c[["jagsresult"]]$BUGSoutput$summary
split.df <- as.data.frame(split.df[grepl("^d\\[[0-9]+\\]",
                                         rownames(split.df)), c(3, 5, 7)])
split.df$treatment <- nma.c[["trt.labs"]]
split.df$agent <- sapply(nma.c[["trt.labs"]], function(x) strsplit(x,
                                                                   split = "_", fixed = TRUE)[[1]][1])
split.df$dose <- as.numeric(sapply(nma.c[["trt.labs"]],
                                   function(x) strsplit(x, split = "_", fixed = TRUE)[[1]][2]))
split.df$agent <- factor(split.df$agent, levels = network$agents)
split.df <- split.df[split.df$agent %in% g$data$agent, ]
lE0 <- rescale.link(E0, direction = "link", link = "logit")
split.df[, c(1:3)] <- split.df[, c(1:3)] + lE0
for (i in 1:3) {
  split.df[, i] <- rescale.link(split.df[, i], direction = "natural",
                                link = "logit")
}
split.df$dose[split.df$treatment=="Brodalumab_105"] <- 100  # CHECK THIS!!!!!!!
split.df$Model <- NA


# Plot predicted responses for all models (Figure 5)
g <- ggplot(plot.df, aes(x=dose, y=`50%`, ymin=`2.5%`, ymax=`97.5%`, group=Model, linetype=Model)) +
  geom_ribbon(aes(fill=Model), alpha=0.2) +
  geom_line() +
  facet_wrap(~agent, scales="free_x") +
  geom_point(data = split.df, ggplot2::aes(x = dose, y = `50%`), size=1.5) +
  geom_linerange(data = split.df, aes(x = dose, ymin = `2.5%`, ymax = `97.5%`), linetype="solid", size=0.7) +
  xlab("Dose mg/week") + ylab("Proportion of patients with >=75% PASI improvement") +
  theme_bw()

plot(g)





###########################################################################
################## Run Unrelated Mean Effect (UME) Models #################
###########################################################################

################## Run NMA UME Models #################

# Common effects treatment-level "split" UME model
nma.ume.c <- nma.run(network, n.iter=15000, pd="pd.kl", n.chain=4, UME = TRUE,
                     parameters.to.save=c("d", "totresdev", "resdev"))

# Random effects treatment-level "split" UME model
nma.ume.r <- nma.run(network, method="random", n.iter=50000, n.burnin=30000,
                     pd="pd.kl", n.chain=4, UME=TRUE,
                     parameters.to.save=c("d", "totresdev", "resdev", "sd"))

# Common effects agent-level "lumped" UME model
nma.agent.ume.c <- nma.run(netagent, n.iter=50000, n.burnin=30000, pd="pd.kl",
                           n.chain=4, UME=TRUE,
                           parameters.to.save=c("d", "totresdev", "resdev"))

# Random effects agent-level "lumped" UME model
nma.agent.ume.r <- nma.run(netagent, method="random", n.iter=50000, n.burnin=30000,
                           pd="pd.kl", n.chain=4, UME=TRUE,
                           parameters.to.save=c("d", "totresdev", "resdev", "sd"))



############ Run MBNMA UME Models ##############

# Common effects exponential MBNMA UME model
exp.ume.c <- mbnma.run(network, fun=dexp(), pd="pd.kl", n.iter=50000, n.burnin=30000,
                               n.chain=4, UME=TRUE,
                               parameters.to.save=c("rate", "totresdev", "resdev"))

# Random effects exponential MBNMA UME model
exp.ume.r <- mbnma.run(network, fun=dexp(), method="random", pd="pd.kl", n.iter=50000,
                               n.burnin=30000, n.chain=4, UME=TRUE,
                               parameters.to.save=c("rate", "totresdev", "resdev", "sd"))

# Common effects Emax MBNMA UME model
emax.ume.c <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=100000, n.burnin=50000,
                         n.chain=4, UME=TRUE,
                         parameters.to.save=c("emax", "ed50", "totresdev", "resdev"))

# Random effects Emax MBNMA UME model (non-converging)
emax.ume.r <- mbnma.run(network, fun=demax(), pd="pd.kl", n.iter=100000, n.burnin=50000,
                         n.chain=4, method="random", UME=TRUE,
                         parameters.to.save=c("emax", "ed50", "totresdev", "resdev", "sd"))



###########################################################
########### Table of Model Results (Table 1) ##############
###########################################################

neatcri <- function(sumrow, digits=2) {
  str <- paste0(
    round(sumrow[5], digits), " (",
    round(sumrow[3], digits), ", ",
    round(sumrow[7], digits), ")"
  )
  return(str)
}


out.df <- data.frame(
  "Model"=c(rep("NMA treatment-level",2),
            rep("NMA agent-level", 2),
            rep("MBNMA exponential",2),
            rep("MBNMA Emax",2),
            rep("UME treatment-level",2)
            ),
  "TreatmentEffect"=c(
    rep(c("Common", "Random"), 5)
  ),
  "ResidualDeviance"=c(
    nma.c$jagsresult$BUGSoutput$mean$totresdev,
    nma.r$jagsresult$BUGSoutput$mean$totresdev,
    nma.agent.c$jagsresult$BUGSoutput$mean$totresdev,
    nma.agent.r$jagsresult$BUGSoutput$mean$totresdev,
    exp.c$BUGSoutput$mean$totresdev,
    exp.r$BUGSoutput$mean$totresdev,
    emax.c$BUGSoutput$mean$totresdev,
    emax.r$BUGSoutput$mean$totresdev,
    nma.ume.c$jagsresult$BUGSoutput$mean$totresdev,
    nma.ume.r$jagsresult$BUGSoutput$mean$totresdev
  ),
  "pD"=c(
    nma.c$jagsresult$BUGSoutput$pD,
    nma.r$jagsresult$BUGSoutput$pD,
    nma.agent.c$jagsresult$BUGSoutput$pD,
    nma.agent.r$jagsresult$BUGSoutput$pD,
    exp.c$BUGSoutput$pD,
    exp.r$BUGSoutput$pD,
    emax.c$BUGSoutput$pD,
    emax.r$BUGSoutput$pD,
    nma.ume.c$jagsresult$BUGSoutput$pD,
    nma.ume.r$jagsresult$BUGSoutput$pD
  ),
  "DIC"=c(
    nma.c$jagsresult$BUGSoutput$DIC,
    nma.r$jagsresult$BUGSoutput$DIC,
    nma.agent.c$jagsresult$BUGSoutput$DIC,
    nma.agent.r$jagsresult$BUGSoutput$DIC,
    exp.c$BUGSoutput$DIC,
    exp.r$BUGSoutput$DIC,
    emax.c$BUGSoutput$DIC,
    emax.r$BUGSoutput$DIC,
    nma.ume.c$jagsresult$BUGSoutput$DIC,
    nma.ume.r$jagsresult$BUGSoutput$DIC
  ),
  "SD"=c(
    "-", neatcri(nma.r$jagsresult$BUGSoutput$summary["sd",]),
    "-", neatcri(nma.agent.r$jagsresult$BUGSoutput$summary["sd",]),
    "-", neatcri(exp.r$BUGSoutput$summary["sd",]),
    "-", neatcri(emax.r$BUGSoutput$summary["sd",]),
    "-", neatcri(nma.ume.r$jagsresult$BUGSoutput$summary["sd",])
  ),
  "Selected"=c(
    c("Yes", "No", "No", "Yes", "No", "Yes", "Yes", "No", "-", "-")
  )
)

# Print Table 1
knitr::kable(out.df, digits = 1,
                  col.names = c("Model", "Treatment effect model",
                                "Posterior mean residual deviance", "pD",
                                "DIC", "Between-study SD (95%CrI)", "Selected (based on DIC)"))


############################################################
##### Global Assessment of Consistency (dev-dev plots) #####
############################################################

# Get deviance contributions from NMA models
nmadevs <- function(nma, ume, modnam="NMA treatment-level") {
  dev.nma <- vector()
  dev.ume <- vector()
  study <- vector()
  arm <- vector()

  for (i in 1:nma$jagsresult$model$data()$NS) {
    narm <- nma$jagsresult$model$data()$narm[i]
    dev.nma <- append(dev.nma, nma$jagsresult$BUGSoutput$mean$resdev[i,1:narm])
    dev.ume <- append(dev.ume, ume$jagsresult$BUGSoutput$mean$resdev[i,1:narm])
    study <- append(study, rep(i, narm))
    arm <- append(arm, 1:narm)
  }

  trt.df <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.ume, "study"=study, "arm"=arm, "model"=modnam)
  return(trt.df)
}



# Get deviance contributions from MBNMA models
mbnmadevs <- function(mbnma, ume, nma.ume, modnam="NMA treatment-level") {
  dev.nma <- vector()
  dev.ume <- vector()
  dev.nma.ume <- vector()
  study <- vector()
  arm <- vector()
  for (i in 1:mbnma$model$data()$NS) {
    narm <- mbnma$model$data()$narm[i]
    dev.nma <- append(dev.nma, mbnma$BUGSoutput$median$resdev[i,1:mbnma$model$data()$narm[i]])
    dev.ume <- append(dev.ume, ume$BUGSoutput$median$resdev[i,1:ume$model$data()$narm[i]])
    dev.nma.ume <- append(dev.nma.ume, nma.ume$jagsresult$BUGSoutput$median$resdev[i,1:nma.ume$jagsresult$model$data()$narm[i]])

    study <- append(study, rep(i, narm))
    arm <- append(arm, 1:narm)
  }

  trt.df <- data.frame("nma"=dev.nma, "ume"=dev.ume, "nma.ume"=dev.nma.ume, "study"=study, "arm"=arm,
                       "model"= modnam)

  return(trt.df)
}




##### Generate data frame of deviance contributions from all models #######

# Common effects treatment-level NMA
trt.df <- nmadevs(nma.c, nma.ume.c, modnam="NMA treatment-level (common effects)")

# Random effects treatment-level NMA
trt.df <- rbind(trt.df, nmadevs(nma.r, nma.ume.r, modnam="NMA treatment-level (random effects)"))

# Common effects agent-level NMA
trt.df <- rbind(trt.df, nmadevs(nma.agent.c, nma.agent.ume.c, modnam="NMA agent-level (common effects)"))

# Random effects agent-level NMA
trt.df <- rbind(trt.df, nmadevs(nma.agent.r, nma.agent.ume.r, modnam="NMA agent-level (random effects)"))

# Common effects exponential MBNMA
trt.df <- rbind(trt.df, mbnmadevs(exp.c, exp.ume.c, nma.ume = nma.ume.c, modnam="Exponential MBNMA (common effects)"))

# Random effects exponential MBNMA
trt.df <- rbind(trt.df, mbnmadevs(exp.r, exp.ume.r, nma.ume = nma.ume.r, modnam="Exponential MBNMA (random effects)"))

# Common effects Emax MBNMA
trt.df <- rbind(trt.df, mbnmadevs(emax.c, emax.ume.c, nma.ume = nma.ume.c, modnam="Emax MBNMA (common effects)"))

# Random effects Emax MBNMA
trt.df <- rbind(trt.df, mbnmadevs(emax.r, emax.ume.r, nma.ume = nma.ume.r, modnam="Emax MBNMA (random effects)"))


# Generate model labels
dev.df <- trt.df %>% mutate(model=factor(model))

# Define whether "reference" model has common or random treatment effects
dev.df$split <- dev.df$nma[dev.df$model=="NMA treatment-level (common effects)"]
dev.df$split[grepl("random", dev.df$model)] <- dev.df$nma[dev.df$model=="NMA treatment-level (random effects)"]




############ Common Effect models (Figure 6) ##############

# Generate indicator for points to highlight/star
dev.df$incon <- 0
dev.df$incon[dev.df$study==16 & dev.df$arm==1] <- 1
dev.df$incon[dev.df$study==15 & dev.df$arm==1] <- 1
dev.df$incon <- factor(dev.df$incon)

mods <- c("Emax MBNMA (common effects)",
          "Exponential MBNMA (common effects)",
          "NMA agent-level (common effects)",
          "NMA treatment-level (common effects)")

# Generate panel a)
g1 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=split, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: NMA treatment-level model") +
  xlim(c(0,15)) + ylim(c(0,15)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_bw() +
  theme(legend.position = "none")


# Generate panel b)
g2 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=nma.ume, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: UME treatment-level model") +
  xlim(c(0,15)) + ylim(c(0,15)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_bw() +
  theme(legend.position = "none")


# Plot common effect models (Figure 6)
cowplot::plot_grid(g1, g2, nrow=2, labels=c("a)", "b)"), label_size = 20,
                   label_fontfamily = "serif",
                   label_fontface = "plain")





############ Random Effect models (Figure S2) ##############

# Generate indicator for points to highlight/star
dev.df$incon <- 0
dev.df$incon[dev.df$study==16 & dev.df$arm==1] <- 1
dev.df$incon[dev.df$study==15 & dev.df$arm==1] <- 1
dev.df$incon <- factor(dev.df$incon)

mods <- c("Emax MBNMA (common effects)",
          "Exponential MBNMA (random effects)",
          "NMA agent-level (random effects)",
          "NMA treatment-level (common effects)")

# Generate panel a)
g1 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=split, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: NMA treatment-level model") +
  xlim(c(0,15)) + ylim(c(0,15)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_bw() +
  theme(legend.position = "none")


# Generate panel b)
g2 <- ggplot(dev.df[dev.df$model %in% mods,],
             aes(y=nma.ume, x=nma, color=incon, shape=incon)) +
  geom_point() +
  geom_abline(color="red") +
  facet_wrap(~ factor(model), scales = "fixed", ncol=2, nrow=2) +
  xlab("Deviance: consistency model") + ylab("Deviance: UME treatment-level model") +
  xlim(c(0,15)) + ylim(c(0,15)) +
  scale_color_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(16,8)) +
  theme_bw() +
  theme(legend.position = "none")


# Plot common effect models (Figure S2)
cowplot::plot_grid(g1, g2, nrow=2, labels=c("a)", "b)"), label_size = 20,
                   label_fontfamily = "serif",
                   label_fontface = "plain")




######################################
##### Node-Splitting  #####
######################################

###### Run selected node-split models #######

# Identify potential node-split comparisons
comps <- inconsistency.loops(network$data.ab, incldr = TRUE)


# Run node-split NMA models (treatment-level) (common effects selected over random effects)
split.nma <- nma.nodesplit(network, comparisons = comps[c(1,3:5),])

# Node-split exponential MBNMA (random effects selected over common effects)
split.exp <- mbnma.nodesplit(network, fun=dexp(), method="random", n.iter=50000, n.burnin=30000, n.chain=4,
                               comparisons = comps)

# Node-split Emax MBNMA (common effects selected over random effects)
split.emax <- mbnma.nodesplit(network, fun=demax(), n.iter=200000, n.burnin=70000, n.chain=4,
                              comparisons = comps[1:5,])




########### Plot forest plot of nodesplit results (Figure 7) ###########

library(forestplot)

nma.ns <- summary(split.nma)
nma.ns$model <- "Treatment-level NMA"

exp.ns <- summary(split.exp)
exp.ns$model <- "Exponential MBNMA"

emax.ns <- summary(split.emax)
emax.ns$model <- "Emax MBNMA"

# Comparisons with indirect evidence via consistency
cons.ns.1 <- rbind(nma.ns,
                   emax.ns[emax.ns$Comparison %in% c(nma.ns$Comparison),],
                   exp.ns[exp.ns$Comparison %in% c(nma.ns$Comparison),]
)

# Comparisons with indirect evidence only via dose-response (>=3 doses)
temp <- emax.ns[!emax.ns$Comparison %in% c(nma.ns$Comparison),]
cons.ns.2 <- rbind(temp,
                   exp.ns[exp.ns$Comparison %in% c(temp$Comparison),]
)

# Comparisons with indirect evidence only via dose-response (>=2 doses)
cons.ns.3 <- exp.ns[!exp.ns$Comparison %in% c(emax.ns$Comparison),]



# Create data frame of node-split results from all models
cons.ns <- rbind(cons.ns.1, cons.ns.2, cons.ns.3)

# Define evidence contribution
cons.ns$Evidence[cons.ns$Evidence %in% c("MBNMA", "NMA")] <- "Overall"

cons.ns <- cons.ns %>% group_by(Comparison) %>% mutate(count=n())
cons.ns <- arrange(cons.ns, desc(count), Comparison, model, Evidence)

names(cons.ns)[4:5] <- c("l95", "u95")

# Write comparison labels
cons.ns$Comparison <- as.character(factor(cons.ns$Comparison, labels=c(
  expression("Ix 40 vs\nIx 20"),
  expression("Us 3.75 vs\nEt 100"),
  expression("Us 5.62 vs\nSe 75"),
  expression("Us 7.5 vs\nEt 100"),

  expression("Us 5.62 vs\nIx 30"),

  expression("Se 37.5 vs\nPl 0"),
  expression("Ti 8.3 vs\nPl 0"),
  expression("Ti 16.7 vs\nPl 0"),
  expression("Br 70 vs\nPl 0"),
  expression("Br 105 vs\nPl 0"),
  expression("Et 50 vs\nPl 0"),
  expression("In 34.5 vs\nPl 0"),
  expression("In _57.6 vs\nPl 0")
)))


# Add column for p-values and format correctly
cons.ns$p.value <- format(cons.ns$p.value)
cons.ns$p.value[cons.ns$p.value=="0.000"] <- "<0.001"


# Remove surplus labels
cons.ns <- cons.ns %>% group_by(Comparison, model) %>% mutate(count=seq(n()))
cons.ns$model[cons.ns$count!=2] <- NA
cons.ns$p.value[cons.ns$count!=2] <- NA

cons.ns <- cons.ns %>% group_by(Comparison) %>% mutate(count=seq(n()))
cons.ns$Comparison[cons.ns$count>1] <- NA


# Add blank between models for plot spacing
index <- which(cons.ns$Evidence=="Overall")

temp <- cons.ns[1:(index[1]),]
for (i in 1:(length(index)-1)) {
  temp[nrow(temp)+1,] <- NA
  temp <- rbind(temp, cons.ns[(index[i]+1):index[i+1],])
}
cons.ns <- temp


# Add blank between Comparisons for plot spacing
index <- which(!is.na(cons.ns$Comparison))
index <- c(index, nrow(cons.ns)+1)

temp <- cons.ns[0,]
index <- index[]
for (i in 1:(length(index)-1)) {
  temp[nrow(temp)+1,] <- NA
  temp <- rbind(temp, cons.ns[(index[i]):(index[i+1]-1),])
}
temp[nrow(temp)+1,] <- NA
cons.ns <- temp


# Add dividing lines between each comparison
index <- which(!is.na(cons.ns$Comparison))[-1]
sumind <- rep(FALSE, nrow(cons.ns))
sumind[index-2] <- TRUE

cex <- 0.6 # Set text size

# Plot forest plot of node-split results (Figure 7)
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





