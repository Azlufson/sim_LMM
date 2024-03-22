library(tidyr)
library(dplyr)
library(lme4)
library(lmerTest)


##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#n.subj und n.obs müssen gerade sein
sim_data_int <- function(n.subj = 10, n.obs = 6, b0 = 10, beta_obs = 0, beta_cond = 5, sd.int_subj = 6, sd_eps = 2) {
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
#   data <- sim_data_int(n.subj = n.subj, n.obs = n.obs, b0 = 10, beta_obs = beta_obs, beta_cond = 5, sd.int_subj = 5, sd_eps = 1)
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

#full and null model (LRT):
m.full <- y ~ obs + cond + (1|subj)
m.null <- y ~ cond + (1|subj)

#model
model <- y ~ obs + cond + (1|subj)

#Parameter für Simulationen
nsim <- 2500
if(nsim < 2) nsim <- 2 #sonst funktioniert Skript nicht
beta_obs <- 0 #auf diesen fixed effect wird jeweils getestet
n.subj <- c(4, 6, 10, 16)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(n.subj, n.obs)
colnames(grid) <- c("n.subj", "n.obs")

#Seed
set.seed(1996)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(n.subj = x[1], n.obs = x[2], beta_obs = beta_obs), m.full, m.null))))
data_LRT.REML_long <- cbind(grid, data_LRT.REML)
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, (ncol(grid)+1):ncol(data_LRT.REML_long))

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
data_LRT.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(n.subj = x[1], n.obs = x[2], beta_obs = beta_obs), m.full, m.null, REML = FALSE))))
data_LRT.ML_long <- cbind(grid, data_LRT.ML)
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, (ncol(grid)+1):ncol(data_LRT.ML_long))

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
data_TasZ.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_obs = beta_obs), model))))
data_TasZ.REML_long <- cbind(grid, data_TasZ.REML)
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, (ncol(grid)+1):ncol(data_TasZ.REML_long))

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
data_TasZ.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_obs = beta_obs), model, REML = FALSE))))
data_TasZ.ML_long <- cbind(grid, data_TasZ.ML)
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, (ncol(grid)+1):ncol(data_TasZ.ML_long))

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
data_SW.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_obs = beta_obs), model, REML = REML, ddf = ddf))))
data_SW.REML_long <- cbind(grid, data_SW.REML)
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, (ncol(grid)+1):ncol(data_SW.REML_long))

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
data_KR.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_obs = beta_obs), model, REML = TRUE, ddf = ddf))))
data_KR.REML_long <- cbind(grid, data_KR.REML)
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, (ncol(grid)+1):ncol(data_KR.REML_long))

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
data_SW.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(n.subj = x[1], n.obs = x[2], beta_obs = beta_obs), model, REML = TRUE, ddf = ddf))))
data_SW.ML_long <- cbind(grid, data_SW.ML)
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, (ncol(grid)+1):ncol(data_SW.ML_long))

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

### Grafiken der Ergebnisse
data_n <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML)
data_n$n.obs <- as.factor(data_n$n.obs)
data_n$n.subj <- as.factor(data_n$n.subj)
data_n$REML <- factor(data_n$REML, labels = c("ML", "REML"))
data_n$method <- factor(data_n$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger"))

save(data_n, file = paste0("/data/hermann/data_n_", nsim, ".RData"))
rm(list = ls())


#############################################################################################################################
#Missing Values
#############################################################################################################################

##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#n.subj und n.obs müssen gerade sein
sim_data_int <- function(n.subj = 10, n.obs = 6, b0 = 10, beta_obs = 0, beta_cond = 5, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs)
  obs <- rep(rep(c(0,1), each = n.obs/2), n.subj)
  cond <- rep(c(0,1), n.subj*n.obs/2)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs)
  y <- b0 + beta_obs * obs + beta_cond * cond + subj_int + rnorm(n.obs*n.subj, 0, sd_eps)
  return(data.frame(subj, obs, cond, y))
}

##missing values erzeugen (MCAR)
missing.val <- function(data, p) {
  return(data[-sample(1:nrow(data), size = p * nrow(data), replace = FALSE),])
}

test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), df = 1, lower = FALSE))
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
nsim <- 2500
beta_obs <- 0 #auf diesen fixed effect wird jeweils getestet
p.missing <- c(.1, .3, .5)

set.seed(1996)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(sapply(p.missing, function(x) replicate(nsim, test_lrtstat(missing.val(sim_data_int(n.subj = 100), p = x), m.full, m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(p.missing, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 2:ncol(data_LRT.REML_long))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.LRT.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(sapply(p.missing, function(x) replicate(nsim, test_lrtstat(missing.val(sim_data_int(), p = x), m.full, m.null, REML = FALSE))))
colnames(data_LRT.ML) <- 1:nsim
data_LRT.ML_long <- as.data.frame(cbind(p.missing, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 2:ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.LRT.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(sapply(p.missing, function(x) replicate(nsim, test_TasZ.fixed(missing.val(sim_data_int(), p = x), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(p.missing, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 2:ncol(data_TasZ.REML_long))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(sapply(p.missing, function(x) replicate(nsim, test_TasZ.fixed(missing.val(sim_data_int(), p = x), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(p.missing, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 2:ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(sapply(p.missing, function(x) replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(p.missing, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 2:ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.SW.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(sapply(p.missing, function(x) replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = TRUE, ddf = ddf))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(p.missing, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 2:ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.KR.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(sapply(p.missing, function(x) replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = TRUE, ddf = ddf))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(p.missing, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 2:ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.SW.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!

### Grafiken der Ergebnisse
data_missing <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML)
data_missing$p.missing <- as.factor(data_missing$p.missing)
data_missing$REML <- factor(data_missing$REML, labels = c("ML", "REML"))
data_missing$method <- factor(data_missing$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger"))

save(data_missing, file = paste0("/data/hermann/data_missing_", nsim, ".RData"))
rm(list = ls())

#############################################################################################################################
#Unbalanced Data
#############################################################################################################################

##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + beta_group * group + (1|cond)
#n.subj und n.obs müssen gerade sein
sim_data_int_unb <- function(n.group1 = 10, n.group2 = 10, n.cond = 4, b0 = 10, beta_group = 0, sd.int_cond = 6, sd_eps = 1) {
  group <- c(rep(1, n.group1 * n.cond), rep(2, n.group2 * n.cond))
  cond <- rep(1:n.cond, length(group)/n.cond)
  cond_int <- rep(rnorm(n.cond, 0, sd.int_cond), length(group)/4)
  y <- b0 + beta_group * group + cond_int + rnorm(length(group), 0, sd_eps)
  return(data.frame(group = as.factor(group), cond = as.factor(cond), y))
}

test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), df = 1, lower = FALSE))
  #return(as.numeric(2 * (logLik(full) - logLik(null))))
}

##t-as-z
test_TasZ.fixed <- function(data, m.full, REML = TRUE) {
  return(summary(lmer(model, data = data, REML = REML))$coefficients[2,4])
}

##KR, SW
#ddf ... Art der Approximation
test_approx.fixed <- function(data, model, REML = TRUE, ddf) {
  return(anova(lmer(model, data = data, REML = REML), ddf = ddf)$`Pr(>F)`[1])
}

#full and null model (LRT):
m.full <- y ~ group + (1|cond)
m.null <- y ~ (1|cond)

#model
model <- y ~ group + (1|cond)

#Parameter für Simulationen
nsim <- 2500
beta_obs <- 0 #auf diesen fixed effect wird jeweils getestet
n.group1 <-  10
n.group2 <- n.group1 * c(1, 2, 3, 4)
grid <- expand.grid(n.group1, n.group2)
colnames(grid) <- c("n.group1", "n.group2")

set.seed(1996)

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(data = sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full = m.full, m.null = m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(grid, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 3:ncol(data_LRT.REML_long))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.LRT.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full, m.null, REML = FALSE))))
colnames(data_LRT.ML) <- 1:nsim
data_LRT.ML_long <- as.data.frame(cbind(grid, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 3:ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.LRT.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(grid, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 3:ncol(data_TasZ.REML_long))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(grid, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 3:ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(grid, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 3:ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.SW.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(grid, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 3:ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.KR.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(grid, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 3:ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.SW.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!


### Grafiken der Ergebnisse
data_unbalanced <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML)
data_unbalanced$n.group1 <- as.factor(data_unbalanced$n.group1)
data_unbalanced$n.group2 <- as.factor(data_unbalanced$n.group2)
data_unbalanced$REML <- factor(data_unbalanced$REML, labels = c("ML", "REML"))
data_unbalanced$method <- factor(data_unbalanced$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger"))

save(data_unbalanced, file = paste0("/data/hermann/data_unbalanced_", nsim, ".RData"))
rm(list = ls())

#############################################################################################################################
#Effect Size
#############################################################################################################################

##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#n.subj und n.obs müssen gerade sein
sim_data_int <- function(n.subj = 16, n.obs = 16, b0 = 10, beta_obs = 0, beta_cond = 5, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs)
  obs <- rep(rep(c(0,1), each = n.obs/2), n.subj)
  cond <- rep(c(0,1), n.subj*n.obs/2)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs)
  y <- b0 + beta_obs * obs + beta_cond * cond + subj_int + rnorm(n.obs*n.subj, 0, sd_eps)
  return(data.frame(subj, obs, cond, y))
}

test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), df = 1, lower = FALSE))
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


#full and null model (LRT):
m.full <- y ~ obs + cond + (1|subj)
m.null <- y ~ cond + (1|subj)

#model
model <- y ~ obs + cond + (1|subj)

#Parameter für Simulationen
nsim <- 2500
beta_obs <- 0 #auf diesen fixed effect wird jeweils getestet
ES <- seq(0, 1.4, .2)

set.seed(1996)


###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(sapply(ES, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_obs = x), m.full, m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(ES, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 2:ncol(data_LRT.REML_long))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(p.LRT.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(sapply(ES, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_obs = x), m.full, m.null, REML = FALSE))))
colnames(data_LRT.ML) <- 1:nsim
data_LRT.ML_long <- as.data.frame(cbind(ES, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 2:ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(p.LRT.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(sapply(ES, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_obs = x), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(ES, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 2:ncol(data_TasZ.REML_long))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(sapply(ES, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_obs = x), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(ES, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 2:ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(sapply(ES, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_obs = x), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(ES, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 2:ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(p.SW.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(sapply(ES, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_obs = x), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(ES, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 2:ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(p.KR.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(sapply(ES, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_obs = x), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(ES, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 2:ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(ES) %>% 
  summarize(k = sum(p.SW.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!

### Grafiken der Ergebnisse
data_ES <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML)
data_ES$ES <- as.factor(data_ES$ES)
data_ES$REML <- factor(data_ES$REML, labels = c("ML", "REML"))
data_ES$method <- factor(data_ES$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger"))

save(data_ES, file = paste0("/data/hermann/data_ES_", nsim, ".RData"))
rm(list = ls())


#############################################################################################################################
#Effect Size by Sample Size
#############################################################################################################################


##Datengeneration
#einfaches Modell nur mit random intercept
# y = b0 + b1*obs + b2*cond + (1|subj) + epsilon
#n.subj und n.obs müssen gerade sein
sim_data_int <- function(n.subj = 16, n.obs = 16, b0 = 10, beta_obs = 0, beta_cond = 5, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs)
  obs <- rep(rep(c(0,1), each = n.obs/2), n.subj)
  cond <- rep(c(0,1), n.subj*n.obs/2)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs)
  y <- b0 + beta_obs * obs + beta_cond * cond + subj_int + rnorm(n.obs*n.subj, 0, sd_eps)
  return(data.frame(subj, obs, cond, y))
}

test_lrtstat <- function(data, m.full, m.null, REML = TRUE) {
  full <- lmer(m.full, data = data, REML = REML)
  null <- lmer(m.null, data = data, REML = REML)
  return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), df = 1, lower = FALSE))
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


#full and null model (LRT):
m.full <- y ~ obs + cond + (1|subj)
m.null <- y ~ cond + (1|subj)

#model
model <- y ~ obs + cond + (1|subj)

#Parameter für Simulationen
nsim <- 2500
beta_obs <- 0 #auf diesen fixed effect wird jeweils getestet
ES <- seq(0, 1.4, .2)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(ES, n.obs)
colnames(grid) <- c("ES", "n.obs")

set.seed(1996)


###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_obs = x[1], n.obs = x[2]), m.full, m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(grid, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 3:ncol(data_LRT.REML_long))

##Daten für Plot:
p_LRT.REML <- data_LRT.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.LRT.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 1)

##ML
data_LRT.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_obs = x[1], n.obs = x[2]), m.full, m.null, REML = FALSE))))
data_LRT.ML_long <- as.data.frame(cbind(grid, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 3:ncol(data_LRT.ML_long))

p_LRT.ML <- data_LRT.ML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.LRT.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 1)

###t-as-z
##REML
data_TasZ.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_obs = x[1], x[2]), model))))
data_TasZ.REML_long <- as.data.frame(cbind(grid, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 3:ncol(data_TasZ.REML_long))

p_TasZ.REML <- data_TasZ.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.REML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 2)

##ML
data_TasZ.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_obs = x[1], n.obs = x[2]), model, REML = FALSE))))
data_TasZ.ML_long <- as.data.frame(cbind(grid, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 3:ncol(data_TasZ.ML_long))

p_TasZ.ML <- data_TasZ.ML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(abs(p.TasZ.ML) >= 1.96) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 2)

###KR, SW
#anova aus lmertest

##Sattherthwaire, REML
ddf <- "Satterthwaite"
REML <- TRUE
data_SW.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_obs = x[1], n.obs = x[2]), model, REML = REML, ddf = ddf))))
data_SW.REML_long <- as.data.frame(cbind(grid, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 3:ncol(data_SW.REML_long))

p_SW.REML <- data_SW.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.SW.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 3)

##Kenward-Roger, REML
ddf <- "Kenward-Roger"
REML <- TRUE
data_KR.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_obs = x[1], n.obs = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
data_KR.REML_long <- as.data.frame(cbind(grid, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 3:ncol(data_KR.REML_long))

p_KR.REML <- data_KR.REML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.KR.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 1, 
         method = 4)

##Sattherthwaire, ML
ddf <- "Satterthwaite"
REML <- FALSE
data_SW.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_obs = x[1], n.obs = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
data_SW.ML_long <- as.data.frame(cbind(grid, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 3:ncol(data_SW.ML_long))

p_SW.ML <- data_SW.ML_long %>% 
  group_by(ES, n.obs) %>% 
  summarize(k = sum(p.SW.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(ES, n.obs, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 3)

#Kenward-Roger nur für ML möglich!

### Grafiken der Ergebnisse
data_ESn <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML)
data_ESn$ES <- as.factor(data_ESn$ES)
data_ESn$n.obs <- as.factor(data_ESn$n.obs)
data_ESn$REML <- factor(data_ESn$REML, labels = c("ML", "REML"))
data_ESn$method <- factor(data_ESn$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger"))

save(data_ESn, file = paste0("/data/hermann/data_ESn_", nsim, ".RData"))
rm(list = ls())