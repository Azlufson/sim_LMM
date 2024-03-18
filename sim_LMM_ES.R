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

####Effektstärke (power?)


##TODO: KR für REML?, PB für REML?
##      funktionen aus anderem skript importieren?
##      ES für ungerade zahlen
##      future_sapply vs future_replicate testen

library(future.apply)
library(parallel)
library(tidyverse)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(afex)

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
#   data <- sim_data_int(n.subj = ES = n.obs, b0 = 10, beta_cond = beta_cond, beta_cond = 5, sd.int_subj = 5, sd_eps = 1)
#   full <- lmer(y ~ obs + cond + (1|subj), data = data, REML = REML)
#   null <- lmer(y ~ cond + (1|subj), data = data, REML = REML)
#   return(pchisq(as.numeric(2 * (logLik(full) - logLik(null))), lower = FALSE))
# }

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
m.full <- y ~ cond + (1|subj)
m.null <- y ~ (1|subj)

#model
model <- y ~ cond + (1|subj)

#Parameter für Simulationen
nsim <- 2
beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
ES <- seq(0, 1.4, .2)
n.obs <- c(4, 6, 10, 16)
grid <- expand.grid(ES, n.obs)
colnames(grid) <- c("ES", "n.obs")

#future_replicate
plan("multisession", workers = detectCores())

#Parameter für parametric bootstrap
nsim.mixed <- 2 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit nsim.pb)
nsim.pb <- 2

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_cond = x[1], n.obs = x[2]), m.full, m.null))))
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
data_LRT.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_lrtstat(sim_data_int(beta_cond = x[1], n.obs = x[2]), m.full, m.null, REML = FALSE))))
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
data_TasZ.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_cond = x[1], x[2]), model))))
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
data_TasZ.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_TasZ.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = FALSE))))
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
data_SW.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = REML, ddf = ddf))))
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
data_KR.REML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
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
data_SW.ML <- t(apply(grid, 1, function(x) replicate(nsim, test_approx.fixed(sim_data_int(beta_cond = x[1], n.obs = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
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
data_ES <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB)
data_ES$ES <- as.factor(data_ES$ES)
data_ES$n.obs <- as.factor(data_ES$n.obs)
data_ES$REML <- factor(data_ES$REML, labels = c("ML", "REML"))
data_ES$method <- factor(data_ES$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

#alle Methoden
ggplot(data_ES, aes(x = ES, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .8) +
  ylim(0, 1)

#nur SW und KR
data_ES %>% 
  filter(method %in% c("Satterthwaite", "Kenward-Roger")) %>% 
  ggplot(aes(x = ES, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .8) +
  ylim(0, 1)

#nur SW
data_ES %>% 
  filter(method %in% c("Satterthwaite")) %>% 
  ggplot(aes(x = ES, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .8) +
  ylim(0, 1)

#nur ML
data_ES %>% 
  filter(REML == "ML") %>% 
  ggplot(aes(x = ES, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .8) +
  ylim(0, 1)

#nur REML
data_ES %>% 
  filter(REML == "REML") %>% 
  ggplot(aes(x = ES, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .8) +
  ylim(0, 1)

#nur t as z
data_ES %>% 
  filter(method == "t-as-z") %>% 
  ggplot(aes(x = ES, y = p, col = REML)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .8) +
  ylim(0, 1)
