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
#lmm <- lmer(model, sim_data_int_unb())
#pvals.fnc(lmm)

####unbalanced designs

#bei unbalaciertem design entstehen viele boundary probleme


##TODO: 
##      fix sim_data_int_unb()
##      KR für REML?, PB für REML?
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
# y = b0 + beta_group * group + (1|cond)
#n.subj und n.obs müssen gerade sein
sim_data_int_unb <- function(n.group1 = 10, n.group2 = 10, n.cond = 4, b0 = 10, beta_group = 0, sd.int_cond = 6, sd_eps = 1) {
  group <- c(rep(0, n.group1 * n.cond), rep(1, n.group2 * n.cond))
  cond <- rep(1:n.cond, length(group)/n.cond)
  cond_int <- rep(rnorm(n.cond, 0, sd.int_cond), length(group)/4)
  y <- b0 + beta_group * group + cond_int + rnorm(length(group), 0, sd_eps)
  return(data.frame(group = as.factor(group), cond = as.factor(cond), y))
}

##Daten unbalanciert machen (nicht mehr benötigt)
#p_rand ... Unbalanciertheit im random effect
#p_fixed ... Unbalanciertheit im fixed effect
# unbalance <- function(data, p_rand = 0, p_fixed) {
#   cond0_ind <- which(data$cond == 0)
#   data$cond[sample(cond0_ind, size = p_rand * length(cond0_ind), replace = FALSE)] <- 1
#   group0_ind <- which(data$group == 0)
#   data$group[sample(group0_ind, size = p_fixed * length(group0_ind), replace = FALSE)] <- 1
#   return(data)
# }

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
m.full <- y ~ group + (1|cond)
m.null <- y ~ (1|cond)

#model
model <- y ~ group + (1|cond)

#Parameter für Simulationen
nsim <- 2
beta_group <- 0 #auf diesen fixed effect wird jeweils getestet
n.group1 <- c(25, 40, 50, 75, 90, 150) 
n.group2 <- c(25, 10, 50, 25, 10, 50)
grid <- cbind(n.group1, n.group2)
#colnames(grid) <- c("n.group1", "n.group2")

plan("multisession", workers = detectCores())

#Parameter für parametric bootstrap
nsim.mixed <- 2 #niedriger, weil pro iteration auch noch gebootstrapped wird (mit nsim.pb)
nsim.pb <- 2

###LRT
##REML (nicht empfohlen)
data_LRT.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_lrtstat(data = sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full = m.full, m.null = m.null))))
colnames(data_LRT.REML) <- 1:nsim
data_LRT.REML_long <- as.data.frame(cbind(grid, data_LRT.REML))
data_LRT.REML_long <- gather(data_LRT.REML_long, sim, p.LRT.REML, 3:ncol(data_LRT.REML_long))

data_LRT.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(prop_LRT.REML = mean(p.LRT.REML <= .05))

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
data_LRT.ML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_lrtstat(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full, m.null, REML = FALSE))))
colnames(data_LRT.ML) <- 1:nsim
data_LRT.ML_long <- as.data.frame(cbind(grid, data_LRT.ML))
data_LRT.ML_long <- gather(data_LRT.ML_long, sim, p.LRT.ML, 3:ncol(data_LRT.ML_long))

data_LRT.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(prop_LRT.ML = mean(p.LRT.ML <= .05))

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
data_TasZ.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_TasZ.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model))))
colnames(data_TasZ.REML) <- 1:nsim
data_TasZ.REML_long <- as.data.frame(cbind(grid, data_TasZ.REML))
data_TasZ.REML_long <- gather(data_TasZ.REML_long, sim, p.TasZ.REML, 3:ncol(data_TasZ.REML_long))

data_TasZ.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(prop_TasZ.REML = mean(abs(p.TasZ.REML) >= 1.96))

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
data_TasZ.ML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_TasZ.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = FALSE))))
colnames(data_TasZ.ML) <- 1:nsim
data_TasZ.ML_long <- as.data.frame(cbind(grid, data_TasZ.ML))
data_TasZ.ML_long <- gather(data_TasZ.ML_long, sim, p.TasZ.ML, 3:ncol(data_TasZ.ML_long))

data_TasZ.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(prop_TasZ.ML = mean(abs(p.TasZ.ML) >= 1.96))

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
data_SW.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = REML, ddf = ddf))))
colnames(data_SW.REML) <- 1:nsim
data_SW.REML_long <- as.data.frame(cbind(grid, data_SW.REML))
data_SW.REML_long <- gather(data_SW.REML_long, sim, p.SW.REML, 3:ncol(data_SW.REML_long))

data_SW.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(prop_SW.REML = mean(p.SW.REML <= .05))

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
data_KR.REML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_KR.REML) <- 1:nsim
data_KR.REML_long <- as.data.frame(cbind(grid, data_KR.REML))
data_KR.REML_long <- gather(data_KR.REML_long, sim, p.KR.REML, 3:ncol(data_KR.REML_long))

data_KR.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(prop_KR.REML = mean(p.KR.REML <= .05))

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
data_SW.ML <- t(apply(grid, 1, function(x) future_replicate(nsim, test_approx.fixed(sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), model, REML = TRUE, ddf = "Satterthwaite"))))
colnames(data_SW.ML) <- 1:nsim
data_SW.ML_long <- as.data.frame(cbind(grid, data_SW.ML))
data_SW.ML_long <- gather(data_SW.ML_long, sim, p.SW.ML, 3:ncol(data_SW.ML_long))

data_SW.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(prop_SW.ML = mean(p.SW.ML <= .05))

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

#future cluster stoppen
future:::ClusterRegistry("stop")

###parametric bootstrap
##ML
REML = FALSE
#Cluster festlegen
n.cores <- parallel::detectCores() - n.cores_no
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#set seed für foreach
registerDoRNG(1996)

data_PB.ML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full = m.full,m.null = m.null, n.bs = n.bs, REML = REML))))
colnames(data_PB.ML) <- 1:nsim.pb
data_PB.ML_long <- as.data.frame(cbind(grid, data_PB.ML))
data_PB.ML_long <- gather(data_PB.ML_long, sim, p.PB.ML, 3:ncol(data_PB.ML_long))

p_PB.ML <- data_PB.ML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.PB.ML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 0,
         method = 5)

##REML
REML = TRUE

#set seed für foreach
registerDoRNG(1996)

data_PB.REML <- t(apply(grid, 1, function(x) replicate(nsim.pb, test_PB(data = sim_data_int_unb(n.group1 = x[1], n.group2 = x[2]), m.full = m.full,m.null = m.null, n.bs = n.bs, REML = REML))))
colnames(data_PB.REML) <- 1:nsim.pb
data_PB.REML_long <- as.data.frame(cbind(grid, data_PB.REML))
data_PB.REML_long <- gather(data_PB.REML_long, sim, p.PB.REML, 3:ncol(data_PB.REML_long))

p_PB.REML <- data_PB.REML_long %>% 
  group_by(n.group1, n.group2) %>% 
  summarize(k = sum(p.PB.REML < .05) + 1.96^2/2,
            n = n() + 1.96^2,
            p = k/n,
            p_l = p - 1.96 * sqrt(p*(1-p)/n),
            p_u = p + 1.96 * sqrt(p*(1-p)/n)) %>% 
  select(n.group1, n.group2, p, p_l, p_u) %>% 
  mutate(REML = 1,
         method = 5)

### Grafiken der Ergebnisse
data_unbalanced <- rbind(p_TasZ.ML, p_TasZ.REML, p_LRT.ML, p_LRT.REML, p_SW.ML, p_SW.REML, p_KR.REML, p_PB.ML, p_PB.REML)
data_unbalanced$n.group1 <- as.factor(data_unbalanced$n.group1)
data_unbalanced$n.group2 <- as.factor(data_unbalanced$n.group2)
data_unbalanced$REML <- factor(data_unbalanced$REML, labels = c("ML", "REML"))
data_unbalanced$method <- factor(data_unbalanced$method, labels = c("LRT", "t-as-z", "Satterthwaite", "Kenward-Roger", "Parametric Bootstrap"))

#alle Methoden
ggplot(data_unbalanced, aes(x = n.group2, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.group1) +
  ylim(0, .12)

#nur SW und KR
data_unbalanced %>% 
  filter(method %in% c("Satterthwaite", "Kenward-Roger")) %>% 
  ggplot(aes(x = n.group2, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.group1, nrow = 1) +
  ylim(0, .1)

#nur SW
data_unbalanced %>% 
  filter(method %in% c("Satterthwaite")) %>% 
  ggplot(aes(x = n.group2, y = p, col = REML, shape = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.group1, nrow = 1) +
  ylim(0, .1)

#nur ML
data_unbalanced %>% 
  filter(REML == "ML") %>% 
  ggplot(aes(x = n.group2, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.group1, nrow = 1) +
  ylim(0, .12)

#nur REML
data_unbalanced %>% 
  filter(REML == "REML") %>% 
  ggplot(aes(x = n.group2, y = p, col = method)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.group1, nrow = 1) +
  ylim(0, .1)

#nur t as z
data_unbalanced %>% 
  filter(method == "t-as-z") %>% 
  ggplot(aes(x = n.group2, y = p, col = REML)) + 
  geom_point(position = position_dodge(.6)) + 
  geom_errorbar(aes(ymin = p_l, ymax = p_u), position = position_dodge(.6), width = .3) +
  geom_hline(yintercept = .05) +
  facet_wrap(~n.group1, nrow = 1) +
  ylim(0, .12)