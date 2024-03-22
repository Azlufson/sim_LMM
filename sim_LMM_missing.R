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

####missing values

##TODO: KR für REML?, PB für REML?
##      funktionen aus anderem skript importieren?
##      ES für ungerade zahlen
##      future_sapply vs future_replicate testen
##      MCAR, MAR, MNAR

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
sim_data_int <- function(n.subj = 10, n.obs = 10, b0 = 10, beta_cond = 0, sd.int_subj = 6, sd_eps = 2) {
  subj <- rep(1:n.subj, each = n.obs * 2)
  cond <- rep(c(0,1), n.subj*n.obs)
  subj_int <- rep(rnorm(n.subj, 0, sd.int_subj), each = n.obs*2)
  y <- 10 + beta_cond * cond + subj_int + rnorm(length(subj), 0, 2)
  return(data.frame(subj, cond, y))
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
m.full <- y ~ cond + (1|subj)
m.null <- y ~ (1|subj)

#model
model <- y ~ cond + (1|subj)

#Parameter für Simulationen
nsim <- 1000
beta_cond <- 0 #auf diesen fixed effect wird jeweils getestet
p.missing <- c(.1, .3, .5)

plan("multisession", workers = detectCores())

#Parameter für parametric bootstrap
nsim.mixed <- 2 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit nsim.pb)
nsim.pb <- 2

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_lrtstat(missing.val(sim_data_int(n.subj = 100), p = x), m.full, m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(p.missing, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 2:ncol(data_LRT.REML_long))

data_LRT.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_LRT.REML = mean(p.LRT.REML <= .05))

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
data_LRT.ML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_lrtstat(missing.val(sim_data_int(), p = x), m.full, m.null, REML = FALSE))))
colnames(data_LRT.ML) <- 1:nsim
data_LRT.ML_long <- as.data.frame(cbind(p.missing, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 2:ncol(data_LRT.ML_long))

data_LRT.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_LRT.ML = mean(p.LRT.ML <= .05))

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
data_TasZ.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_TasZ.fixed(missing.val(sim_data_int(), p = x), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(p.missing, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 2:ncol(data_TasZ.REML_long))

data_TasZ.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_TasZ.REML = mean(abs(p.TasZ.REML) >= 1.96))

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
data_TasZ.ML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_TasZ.fixed(missing.val(sim_data_int(), p = x), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(p.missing, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 2:ncol(data_TasZ.ML_long))

data_TasZ.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_TasZ.ML = mean(abs(p.TasZ.ML) >= 1.96))

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
data_SW.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(p.missing, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 2:ncol(data_SW.REML_long))

data_SW.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_SW.REML = mean(p.SW.REML <= .05))

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
data_KR.REML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = TRUE, ddf = ddf))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(p.missing, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 2:ncol(data_KR.REML_long))

data_KR.REML_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_KR.REML = mean(p.KR.REML <= .05))

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
data_SW.ML <- t(sapply(p.missing, function(x) future_replicate(nsim, test_approx.fixed(missing.val(sim_data_int(), p = x), model, REML = TRUE, ddf = ddf))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(p.missing, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 2:ncol(data_SW.ML_long))

data_SW.ML_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_SW.ML = mean(p.SW.ML <= .05))

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


###parametric bootstrap (nur ML)

#Cluster festlegen (future funktioniert nicht)
nc <- detectCores() # number of cores
cl <- makeCluster(rep("localhost", nc)) # make cluster

data_PB <- t(sapply(p.missing, function(x) replicate(nsim.mixed, test_PB.fixed(model, data = missing.val(sim_data_int(), p = x), nsim.pb = nsim.pb, cl = cl))))
colnames(data_PB) <- 1:nsim.mixed
data_PB_long <- as.data.frame(cbind(p.missing, data_PB))
data_PB_long <- gather(data_PB_long, sim, p.PB, 2:ncol(data_PB_long))

data_PB_long %>% 
  group_by(p.missing) %>% 
  summarize(prop_PB = mean(p.PB <= .05))

p_PB <- data_PB_long %>% 
  group_by(p.missing) %>% 
  summarize(k = sum(p.PB < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(p.missing, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

### Grafiken der Ergebnisse
data_missing <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB)
data_missing$p.missing <- as.factor(data_missing$p.missing)
data_missing$REML <- factor(data_missing$REML, labels = c("ML", "REML"))
data_missing$method <- factor(data_missing$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

#alle Methoden
ggplot(data_missing, aes(x = p.missing, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  ylim(0, .12)

#nur SW und KR
data_missing %>% 
  filter(method %in% c("Satterthwaite", "Kenward-Roger")) %>% 
  ggplot(aes(x = p.missing, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  ylim(0, .1)

#nur SW
data_missing %>% 
  filter(method %in% c("Satterthwaite")) %>% 
  ggplot(aes(x = p.missing, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  ylim(0, .1)

#nur ML
data_missing %>% 
  filter(REML == "ML") %>% 
  ggplot(aes(x = p.missing, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  ylim(0, .12)

#nur REML
data_missing %>% 
  filter(REML == "REML") %>% 
  ggplot(aes(x = p.missing, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  ylim(0, .1)

#nur t as z
data_missing %>% 
  filter(method == "t-as-z") %>% 
  ggplot(aes(x = p,missing, y = p, col = REML)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  ylim(0, .12)