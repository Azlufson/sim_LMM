#Alpha-Fehler, Güte:
#Modell erstellen
#testen: fixed effects, random effects
#Parameter: Stichprobengröße, Stärke des Effekts, Komplexität des Modells, missing values, balanciertes Design
#Schätzmethoden:  ML, REML
#                 Sattherwaite, Kenward-Rogers
#                 MCMC (Baayen et al., 2008)
#                 t as z

#mcmc nicht mehr unterstützt von lme4
#library(languageR)
#lmm <- lmer(model, sim_data_int())
#pvals.fnc(lmm)

####Stichprobengröße

##TODO: KR für REML?, PB für REML?
##      funktionen aus anderem skript importieren?

library(future.apply)
library(parallel)
library(tidyverse)
library(lme4)
library(lmerTest)
library(foreach)
library(doRNG)

##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b2*cond + (1|subj) + epsilon
sim_data_int <- function(n.subj = 10, n.obs = 6, b0 = 10, beta_cond = 0, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs * 2)
  cond <- rep(c(0,1), n.subj*n.obs)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs*2)
  y <- 10 + beta_cond * cond + subj_int + rnorm(length(subj), 0, 2)
  return(data.frame(subj, cond, y))
}

##LRT:
#alt (keine modellspezifikation möglich):
# test_lrtstat.fixed <- function(n.subj = 6, n.obs = 10, beta_cond = 0, REML = TRUE) {
#   data <- sim_data_int(n.subj = n.subj, n.obs = n.obs, b0 = 10, beta_cond = beta_cond, beta_cond = 5, sd.int_subj = 5, sd_eps = 1)
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

##Funktion fürs parametric bootstrapping
#benötigt package foreach
#mixed (bz.w pbmodcomp) parallelisiert nicht
#doRNG notwendig, um richtig seed zu setzen
test_PB <- function(data, m.full, m.null, n.bs, REML = FALSE) {
  nullmod <- lmer(m.null, data = data, REML = REML)
  fullmod <- lmer(m.full, data = data, REML = REML)
  lrstat <- numeric(n.bs)
  lrstat <- foreach(i = 1:n.bs, .combine = "c") %dopar% {
    data$y_sim <- unlist(simulate(nullmod))
    null <- lmer(update(m.null, y_sim ~ .), data = data, REML = REML)
    alt <- lmer(update(m.full, y_sim ~ .), data = data, REML = REML)
    as.numeric(2*(logLik(alt)-logLik(null)))
  }
  return(mean(lrstat > as.numeric(2*(logLik(fullmod)-logLik(nullmod)))))
}

#full and null model (LRT):
m.full <- y ~ cond + (1|subj)
m.null <- y ~ (1|subj)

#model
model <- y ~ cond + (1|subj)

#Parameter für Simulationen
nsim <- 12
beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
n.subj <- c(4, 6, 10, 16)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(n.subj, n.obs)
colnames(grid) <- c("n.subj", "n.obs")

#future_apply
plan("multisession", workers = availableCores(), gc = TRUE)

#Parameter für parametric bootstrap
nsim.pb <- 12 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit n.bs)
n.bs <- 100

#Seed
set.seed(1996)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full, m.null)), future.seed = TRUE))
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
data_LRT.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full, m.null, REML = FALSE)), future.seed = TRUE))
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
data_TasZ.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model)), future.seed = TRUE))
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
data_TasZ.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = FALSE)), future.seed = TRUE))
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
data_SW.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = REML, ddf = ddf)), future.seed = TRUE))
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
data_KR.REML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = TRUE, ddf = ddf)), future.seed = TRUE))
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
data_SW.ML <- t(future_apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), model, REML = TRUE, ddf = ddf)), future.seed = TRUE))
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

#future cluster stoppen
future:::ClusterRegistry("stop")

###parametric bootstrap
##ML
REML = FALSE
#Cluster festlegen
n.cores <- parallel::detectCores()
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#set seed für foreach
registerDoRNG(1996)

data_PB.ML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full = m.full,m.null = m.null, n.bs = n.bs, REML = REML))))
data_PB.ML_long <- cbind(grid, data_PB.ML)
data_PB.ML_long <- gather(data_PB.ML_long, sim, p.PB.ML, (ncol(grid)+1):ncol(data_PB.ML_long))

p_PB.ML <- data_PB.ML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.PB.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

##REML
REML = TRUE

#set seed für foreach
registerDoRNG(123)

data_PB.REML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int(n.subj = x[1], n.obs = x[2], beta_cond = beta_cond), m.full = m.full, m.null = m.null ,n.bs = n.bs, REML = REML))))
parallel::stopCluster(cl = my.cluster)
data_PB.REML_long <- cbind(grid, data_PB.REML)
data_PB.REML_long <- gather(data_PB.REML_long, sim, p.PB.REML, (ncol(grid)+1):ncol(data_PB.REML_long))

p_PB.REML <- data_PB.REML_long %>% 
  group_by(n.subj, n.obs) %>% 
  summarize(k = sum(p.PB.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.obs, n.subj, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 5)

### Grafiken der Ergebnisse
data_n <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB.ML, p_PB.REML)
data_n$n.obs <- as.factor(data_n$n.obs)
data_n$n.subj <- as.factor(data_n$n.subj)
data_n$REML <- factor(data_n$REML, labels = c("ML", "REML"))
data_n$method <- factor(data_n$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

#alle Methoden
ggplot(data_n, aes(x = n.obs, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj) +
  ylim(0, .2)

#nur SW und KR
data_n %>% 
  filter(method %in% c("Satterthwaite", "Kenward-Roger")) %>% 
  ggplot(aes(x = n.obs, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur SW
data_n %>% 
  filter(method %in% c("Satterthwaite")) %>% 
  ggplot(aes(x = n.obs, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur ML
data_n %>% 
  filter(REML == "ML") %>% 
  ggplot(aes(x = n.obs, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .12)

#nur REML
data_n %>% 
  filter(REML == "REML") %>% 
  ggplot(aes(x = n.obs, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .1)

#nur t as z
data_n %>% 
  filter(method == "t-as-z") %>% 
  ggplot(aes(x = n.obs, y = p, col = REML)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.subj, nrow = 1) +
  ylim(0, .12)