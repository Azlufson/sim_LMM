#Alpha-Fehler, Güte:
#Modell erstellen
#testen: fixed effects, random effects
#Parameter: Stichprobengröße, Stärke des Effekts, Komplexität des Modells, missing values, balanciertes Design
#Schätzmethoden:  ML, REML
#                 Sattherwaite, Kenward-Rogers
#                 MCMC (Baayen et al., 2008)
#                 t as z

####Stichprobengröße
#einfaches Modell nur mit random intercept
# y = b0 + b1*obs + b2*cond + (1|subj) + epsilon

##TODO: KR für REML?, PB für REML?
##      funktionen aus anderem skript importieren?
##      implement MCMC

library(future.apply)
library(parallel)
library(tidyverse)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(afex)
library(languageR)

sim_data.n_int <- function(n.subj = 10, n.obs = 6, b0 = 10, beta_obs = 0, beta_cond = 5, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs)
  obs <- rep(rep(c(0,1), each = n.obs/2), n.subj)
  cond <- rep(c(0,1), n.subj*n.obs/2)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs)
  y <- b0 + beta_obs * obs + beta_cond * cond + subj_int + rnorm(n.obs*n.subj, 0, sd_eps)
  return(data.frame(subj, obs, cond, y))
}

##LRT:
#alt (keine modellspezifikation möglich):
# test_lrtstat.fixed <- function(n.subj = 6, n.obs = 10, beta_obs = 0, REML = TRUE) {
#   data <- sim_data.n_int(n.subj = n.subj, n.obs = n.obs, b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1)
#   full <- lmer(y ~ obs + cond + (1|subj), data = data, REML = REML)
#   null <- lmer(y ~ cond + (1|subj), data = data, REML = REML)
#   return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), 1, lower = FALSE))
# }

test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), 1, lower = FALSE))
}

##t-as-z
test_TasZ.fixed <- function(data, m.full, REML = TRUE) {
  return(summary(lmer(model, data = data, REML = REML))$coefficients[2,4])
}

##KR, SW
#ddf ... Art der Approximation
test_approx.fixed <- function(data, model, REML = TRUE, ddf = "Satterthwaite") {
  return(anova(lmer(model, data = data, REML = REML), ddf = ddf)$`Pr(>F)`[1])
}

##Funktion zur Ausgabe des p-Wertes des fixed effects via parametric bootstrap
#mixed auf afex (nutzt pbmodcomp)
#nsim.pb bestimmt anzahl an bootstrap-simulationen von pbmodcomp
#cl erlaubt multicore nutzung (via package parallel)
test_PB.fixed <- function(mode, data, nsim.pb = 1000, cl = NULL) {
  return(suppressMessages(mixed(model, data = data, method = "PB", progress = FALSE, cl = cl, args_test = list(nsim = nsim.pb, cl = cl))$anova_table$`Pr(>PB)`[1]))
}
#suppressMessages: "mixed" will throw a message if numerical variables are not centered on 0, as main effects (of other variables then the numeric one) can be hard to interpret if numerical variables appear in interactions. See Dalal & Zickar (2012).
#kleine Tests haben ergeben, dass es am effizientesten ist, fürs fitten und fürs bootstrap multicore zu nutzen


#full and null model (LRT):
m.full <- y ~ obs + cond + (1|subj)
m.null <- y ~ cond + (1|subj)

#model
model <- y ~ obs + cond + (1|subj)

#Parameter für Simulationen
beta_obs <- 0 #auf diesen fixed effect wird jeweils getestet
n.subj <- c(4, 6, 10, 16)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(n.subj, n.obs)
colnames(grid) <- c("n.subj", "n.obs")

#future_apply
plan("multisession", workers = detectCores())

#Parameter für parametric bootstrap
nsim.mixed <- 100 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit nsim.pb)
nsim.pb <- 500

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1), m.full, m.null)), future.seed = TRUE))
data_LRT.REML_long <- cbind(grid, data_LRT.REML)
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, (ncol(grid)+1):ncol(data_LRT.REML_long))

data_LRT.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_LRT.REML = mean(p.LRT.REML <= .05))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.LRT.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1), m.full, m.null, REML = FALSE)), future.seed = TRUE))
data_LRT.ML_long <- cbind(grid, data_LRT.ML)
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, (ncol(grid)+1):ncol(data_LRT.ML_long))

data_LRT.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_LRT.ML = mean(p.LRT.ML <= .05))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.LRT.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1), model)), future.seed = TRUE))
data_TasZ.REML_long <- cbind(grid, data_TasZ.REML)
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, (ncol(grid)+1):ncol(data_TasZ.REML_long))

data_TasZ.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_TasZ.REML = mean(abs(p.TasZ.REML) >= 1.96))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 2)

##ML
data_TasZ.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1), model, REML = FALSE)), future.seed = TRUE))
data_TasZ.ML_long <- cbind(grid, data_TasZ.ML)
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, (ncol(grid)+1):ncol(data_TasZ.ML_long))

data_TasZ.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_TasZ.ML = mean(abs(p.TasZ.ML) >= 1.96))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1), model, REML = REML, ddf = ddf)), future.seed = TRUE))
data_SW.REML_long <- cbind(grid, data_SW.REML)
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, (ncol(grid)+1):ncol(data_SW.REML_long))

data_SW.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_SW.REML = mean(p.SW.REML <= .05))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.SW.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1), model, REML = TRUE, ddf = "Satterthwaite")), future.seed = TRUE))
data_KR.REML_long <- cbind(grid, data_KR.REML)
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, (ncol(grid)+1):ncol(data_KR.REML_long))

data_KR.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_KR.REML = mean(p.KR.REML <= .05))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.KR.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1), model, REML = TRUE, ddf = "Satterthwaite")), future.seed = TRUE))
data_SW.ML_long <- cbind(grid, data_SW.ML)
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, (ncol(grid)+1):ncol(data_SW.ML_long))

data_SW.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_SW.ML = mean(p.SW.ML <= .05))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.SW.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!


###parametric bootstrap (nur ML)

#Cluster festlegen (future_apply funktioniert nicht)
(nc <- detectCores()) # number of cores
cl <- makeCluster(rep("localhost", nc)) # make cluster

data_PB <- t(apply(grid, 1, function(x) replicate(nsim.mixed, test_PB.fixed(model, data = sim_data.n_int(n.subj = x[1], n.obs = x[2], b0 = 10, beta_obs = beta_obs, beta_cond = 10, sd.int_subj = 10, sd_eps = 1), nsim.pb = nsim.pb, cl = cl))))
data_PB_long <- cbind(grid, data_PB)
data_PB_long <- gather(data_PB_long, sim, p.PB, (ncol(grid)+1):ncol(data_PB_long))

data_PB_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(prop_PB = mean(p.PB <= .05))

p_PB <- data_PB_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.PB < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

### Grafiken der Ergebnisse
data_p <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB)
data_p$n.obs <- as.factor(data_p$n.obs)
data_p$n.subj <- as.factor(data_p$n.subj)
data_p$REML <- factor(data_p$REML, labels = c("ML", "REML"))
data_p$method <- factor(data_p$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

#alle Methoden
ggplot(data_p, aes(x = as.factor(n.obs), y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj) +
  ylim(0, .12)

#nur SW und KR
data_p %>% 
  filter(method %in% c("Satterthwaite", "Kenward-Roger")) %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur SW
data_p %>% 
  filter(method %in% c("Satterthwaite")) %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur ML
data_p %>% 
  filter(REML == "ML") %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .12)

#nur REML
data_p %>% 
  filter(REML == "REML") %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur t as z
data_p %>% 
  filter(method == "t-as-z") %>% 
  ggplot(aes(x = as.factor(n.obs), y = p, col = REML)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .12)