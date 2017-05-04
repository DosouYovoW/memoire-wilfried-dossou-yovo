##########################
# CHARGER LES LIBRAIRIES #
##########################
library(dplyr)# DATA MANIPULATION
library(ggplot2)# GRAPHIQUES
library(compositions)# COMPOSITIONAL ANALYSIS
library(vegan)# ANALYSE EN COMPOSANTES PRINCIPALE
library(VIM) # VISUALISER DONNEES MANQUANTES
library(mice)# IMPUTATION DES DONNEES MANQUANTES
library(nlme)# MODELISATION MULTI-NIVEAUX
source('scripts/ilrNA.R')# ATTRIBUTION DES VALEURS MANQUANTES DANS LES ILR
source('fonctions.R')# GRAPHIQUES MULTI-NIVEAUX


#######################
# CHARGER LES DONNÉES #
#######################
# Tableau des indices climatiques
indices.climatiques = tbl_df(read.csv("donnees/climat_02.csv", sep = ";", dec = ".", 
                                      header = TRUE, na.strings = 'NaN', encoding = 'UTF-8'))

# Tableau des données
data = tbl_df(read.csv("donnees/data_v1.0.csv", sep = ";", dec = ".", header = TRUE, encoding = 'UTF-8'))
data['Localisation'][data['Localisation'] == ""] = NA
data['Localisation'] = droplevels(data['Localisation'])

# Tableau des séries de sol
serie.sol = tbl_df(read.csv("donnees/series_mais_member.csv", sep = ";", dec = ",", header = TRUE, encoding = 'UTF-8'))

# Tableau pour la traduction des termes
translate_col = read.csv(file = "donnees/translate_colnames.csv", header = TRUE, sep = ";", encoding = 'UTF-8')

# Indice UTM des cultivars
UTM_LOC = tbl_df(read.csv("donnees/Localisation.UTM.csv", sep = ";", dec = ".", header = TRUE))[, c(1, 2)]
data = left_join(data, UTM_LOC, by = "Localisation") # joindre les UTM dans le tableau principal

# Enlever les données climatiques imputées arbitrairement, où les coordonnées sont vides
indices.climatiques[data['Coord.Lati'] == "" | data['Coord.Long'] == "", -c(1, 2)] = NA
indices.climatiques['DJCm_.30'] = NULL
indices.climatiques['DJCm'] = NULL

# Créer une variable pour l'amendement organique
levels(data$Dose.FumierOuLisier)
data$amendement.organique = 1
data$amendement.organique[data$Dose.FumierOuLisier == "Aucun"] = 0
data$amendement.organique[data$Dose.FumierOuLisier == "lactosérum"] = 0
data$amendement.organique[data$Dose.FumierOuLisier == ""] = NA

# listes de variables
variable.explicative = c("UTM.FA", "UTM.diff", "Culture.Densite.Plan_plant.ha", "ProprieteSol.pHeauCaCl2",
                         "ProprieteSol.C_.", "ProprieteSol.Argile_.", "ProprieteSol.Sable_.")
### Note: UTM.diff est défini plus loin, et n'est pas nécessaire avant le calcul de data_rs
### "ProprieteSol.MIII.Mn_mg.kg", trop de données manquantes
### "amendement.organique" trop variable, peu informative

variables.dose = c("Dose.P.Quantite_kgP2o5.ha_temoin", "Dose.P.Quantite_kgP2o5.ha_dose1",
                   "Dose.P.Quantite_kgP2o5.ha_dose2", "Dose.P.Quantite_kgP2o5.ha_dose3")

variable.performance = c(
  "Performance.P.Rendement_t.ha_temoin", "Performance.P.Rendement_t.ha_dose1","Performance.P.Rendement_t.ha_dose2", 
  "Performance.P.Rendement_t.ha_dose3", "Performance.P.Rendement_t.ha_dose4",
  "Performance.P.DensiteGrain_g.L_temoin", "Performance.P.DensiteGrain_g.L_dose1", "Performance.P.DensiteGrain_g.L_dose2",
  "Performance.P.DensiteGrain_g.L_dose3", "Performance.P.DensiteGrain_g.L_dose4",
  "Performance.P.HumiditeGrain_._temoin", "Performance.P.HumiditeGrain_._dose1", "Performance.P.HumiditeGrain_._dose2",
  "Performance.P.HumiditeGrain_._dose3", "Performance.P.HumiditeGrain_._dose4",
  "Performance.P.HauteurPlant_cm_temoin", "Performance.P.HauteurPlant_cm_dose1", "Performance.P.HauteurPlant_cm_dose2",
  "Performance.P.HauteurPlant_cm_dose3"
)


#######################
# Indices climatiques #
#######################
# Analyse en composantes principales
## présélection des variables
colnames(indices.climatiques)
any(indices.climatiques$IhBR != indices.climatiques$IhKcBR) # Les chiffres de IhBR sont-ils différents à ceux de IhKcBR
any(indices.climatiques$etpBRcum != indices.climatiques$etpKcBRcum)

## Corrélations entre les variables
### Fonctions pour les pair plots
panel.hist = function(x, ...) {
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h = hist(x, plot = FALSE)
  breaks = h$breaks; nB = length(breaks)
  y = h$counts; y = y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}
panel.cor = function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr = par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r = abs(cor(x, y))
  txt = format(c(r, 0.123456789), digits = digits)[1]
  txt = paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor = 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

### pairs plot
#### Précipitations
pairs(indices.climatiques[c("nbPcp1",
                            "nbPcp2",
                            "nbPcp5",
                            "nbPcp10",
                            "nb3po3j",
                            "nb4po7j",
                            'pcpCum')],
      upper.panel = panel.cor, diag.panel = panel.hist)

pairs(indices.climatiques[c("NbPcp1_.30",
                            "NbPcp2_.30",
                            "NbPcp5_.30",
                            "NbPcp10_.30",
                            "Nb3po3j_.30",
                            "Nb4po7j_.30",
                            'PcpCum_.30')], 
      upper.panel = panel.cor, diag.panel = panel.hist)

#### indices pour l'ACP
pairs(indices.climatiques[complete.cases(indices.climatiques),
                          c('pcpCum',
                            'SDI',
                            'AWDR',
                            'IhBR',
                            'etpBRcum',
                            'tMean',
                            'DJCb',
                            'nb3po3j',
                            'UTMc')],
      upper.panel = panel.cor, diag.panel = panel.hist)

pairs(indices.climatiques[complete.cases(indices.climatiques),
                          c('PcpCum_.30',
                            'SDI_.30',
                            'AWDR_.30',
                            'IhBR_.30',
                            'EtpBRcum_.30',
                            'Tmean_.30',
                            'DJCb_.30',
                            'Nb3po3j_.30',
                            'UTMc_.30')],
      upper.panel = panel.cor, diag.panel = panel.hist)


## ACP
### Saison complète
indices.climatiques %>%
  filter(complete.cases(.))

pca.climat = rda(indices.climatiques[complete.cases(indices.climatiques),
                                     c('pcpCum',
                                       'SDI',
                                       'AWDR',
                                       'IhBR',
                                       'etpBRcum',
                                       'tMean',
                                       'DJCb',
                                       'nb3po3j',
                                       'UTMc')], scale = TRUE)
#summary(pca.climat)
biplot(pca.climat)
text(pca.climat, display = 'species')
data.frame(inertie = 100 * pca.climat$CA$eig / length(pca.climat$CA$eig)) # pourcentage d'intertie des axes

### 0-30 jours
pca.climat_.30 = rda(indices.climatiques[complete.cases(indices.climatiques), 
                                         c('PcpCum_.30',
                                           'SDI_.30',
                                           'AWDR_.30',
                                           'IhBR_.30',
                                           'EtpBRcum_.30',
                                           'Tmean_.30',
                                           'DJCb_.30',
                                           'Nb3po3j_.30',
                                           'UTMc_.30')], scale = TRUE)
biplot(pca.climat_.30)
text(pca.climat_.30, display = 'species')
data.frame(inertie = 100 * pca.climat_.30$CA$eig / length(pca.climat_.30$CA$eig)) # pourcentage d'intertie des axes


selection.climat = c('SDI', 'UTMc', 'pcpCum', 'tMean')
selection.climat_.30 = c('SDI_.30', 'UTMc_.30', 'PcpCum_.30', 'Tmean_.30')
plot(indices.climatiques[selection.climat], upper.panel = panel.cor, diag.panel = panel.hist)

data = cbind(data, indices.climatiques[c(selection.climat, selection.climat_.30)])# FUSION DE DEUX TABLEAUX PAR COLONNE
data["UTM.diff"] = data["UTMc"] - data["UTM.FA"]

##################
# Analyse de sol #
##################
# P / Al
colnames(data)
# Fusion Al AA et Al SEP
data[c("ProprieteSol.MIII.Al.AA_mg.kg",
       "ProprieteSol.MIII.Al.SEP_mg.kg")]
data["ProprieteSol.MIII.Al_mg.kg"] = data["ProprieteSol.MIII.Al.AA_mg.kg"]
isna_Al = is.na(data$ProprieteSol.MIII.Al_mg.kg)
data$ProprieteSol.MIII.Al_mg.kg[isna_Al] = data$ProprieteSol.MIII.Al.SEP_mg.kg[isna_Al]

# Fusion P AA transformé et P SEP
data[c("ProprieteSol.MIII.P.colorimetrie_mg.kg",
       "ProprieteSol.MIII.P.SEP_mg.kg")]
data["ProprieteSol.MIII.P_mg.kg"] = data[c("ProprieteSol.MIII.P.SEP_mg.kg")]

isna_P = is.na(data$ProprieteSol.MIII.P_mg.kg)
data$ProprieteSol.MIII.P_mg.kg[isna_P] = data$ProprieteSol.MIII.P.colorimetrie_mg.kg[isna_P] * 1.09

# valeur de remplissage
Fv_P.Al = 1E6 - (data$ProprieteSol.MIII.Al_mg.kg + data$ProprieteSol.MIII.P_mg.kg)

# SBP Al, P, Fv
sbp_P.Al = matrix(c(-1, 1, 0,
                    1, 1, -1),
                  ncol = 3, byrow = TRUE)
colnames(sbp_P.Al) = c('Al', 'P', 'Fv')
psi_P.Al = gsi.buildilrBase(t(sbp_P.Al))

parts_P.Al = data.frame(
  Al = data$ProprieteSol.MIII.Al_mg.kg,
  P = data$ProprieteSol.MIII.P_mg.kg,
  Fv = Fv_P.Al
)
colnames(parts_P.Al) = c('Al', 'P', 'Fv')
comp_P.Al = acomp(parts_P.Al)
balances_P.AL = ilr(comp_P.Al, V = psi_P.Al)
balances_P.AL = ilrNA(comp = comp_P.Al, sbp = sbp_P.Al, bal = balances_P.AL)
parts_P.Al[1:10,]
balances_P.AL[1:10,]
data$Sol_Al.P = balances_P.AL[, 1]
data$Sol_Fv.AlP = balances_P.AL[, 2]
plot(balances_P.AL, xlab = '[Al | P]', ylab = '[Fv | Al,P]')
plot(x = parts_P.Al$P, y = parts_P.Al$Al)
plot(x = log(parts_P.Al$P), y = log(parts_P.Al$Al))

histogram(data$Sol_Al.P, xlab = '[Al | P]', breaks = 10, col = 'grey80')
histogram(data$Sol_Fv.AlP, xlab = '[Fv | Al,P]', breaks = 10, col = 'grey80')


################
# Série de sol #
################
data = cbind(data, serie.sol[, -c(1, 2)])

soilTypeSBP.3 = matrix(c(-1, 1,-1,
                         -1, 0, 1),
                       byrow = TRUE,# LIRE PAR LIGNE
                       ncol = 3)# NOMBRE DE COLONNE DE MA MATRICE

colnames(soilTypeSBP.3) = c('loam.gley', 'sand.podzol', 'sand.gley')# DONNER UN NOM AU COLONNE
soilTypeComp.3 = acomp(data[c('m.3a', 'm.3b', 'm.3c')])# acomp pour fermer le simplex
soilTypeBal.3 = ilr(soilTypeComp.3, V = gsi.buildilrBase(t(soilTypeSBP.3)))#calcul des ilr
soilTypeBal.3 = ilrNA(comp = soilTypeComp.3, sbp = soilTypeSBP.3, bal = soilTypeBal.3)# MAINTENIR LES NA DE LA MATRICE

data$SerieSol_G.P = soilTypeBal.3[, 1] # ilr1 [gley | podzol]
data$SerieSol_LG.SG = soilTypeBal.3[, 2] # ilr2 gley.[fine | coarse]

histogram(data$SerieSol_G.P, xlab = '[Gley | Podzol]', breaks = 10, col = 'grey80')
histogram(data$SerieSol_LG.SG, xlab = '[Gley loameux | Gley sableux]', breaks = 10, col = 'grey80')

table(data$SerieSol) # compter le nombre d'occurences pour chaque série de sol du tableau

sum(data$SerieSol_LG.SG < 0, na.rm = TRUE) / sum(!is.na(data$SerieSol_LG.SG)) # sum: somme d'un veteur logique. Vrai = 1, Faux = 0.


###########################
# pH et matière organique #
###########################
# Imputation
aggr(data[c('ProprieteSol.pHeauCaCl2', 'ProprieteSol.C_.')])
pH.C_mice = mice(data[c('ProprieteSol.pHeauCaCl2', 'ProprieteSol.C_.')],
                 m = 50, method = 'rf', seed = 6561874)
pH.C_imp = mice::complete(pH.C_mice)
pH.C_imp[is.na(data['ProprieteSol.pHeauCaCl2']) & is.na(data['ProprieteSol.C_.']), ] = NA

## Créer des colonnes avec l'extension _imp pour indiquer qu'elles sont imputées
data['ProprieteSol.pHeauCaCl2_imp'] = pH.C_imp[, 1]
data['ProprieteSol.C_._imp'] = pH.C_imp[, 2]

names(data)[names(data) == 'ProprieteSol.pHeauCaCl2_imp'] = 'pH'

##################
# Texture du sol #
##################
qcTextCentroid = read.csv('donnees/Centroids_texture.csv', sep = ';', dec = '.', header = TRUE)

qcTextCentroid = qcTextCentroid %>%
  filter(Texture %in% levels(data$ProprieteSol.ClasseTexturale)) %>%
  droplevels() # filter to take only the levels present in data$ProprieteSol.ClasseTexturale

colnames(qcTextCentroid)[2] = 'ProprieteSol.ClasseTexturale'

data$ProprieteSol.ClasseTexturale[data$ProprieteSol.ClasseTexturale == 'Loam sablo-graveleux'] = 'Loam sableux'

data = left_join(data, qcTextCentroid, by = "ProprieteSol.ClasseTexturale") # joindre le tableau de droite au tableau de gauche selon la colonne spécifiée dans la colonne by

data$ProprieteSol.Argile_.[is.na(data$ProprieteSol.Argile_.)] = data$ArgileCentroid[is.na(data$ProprieteSol.Argile_.)]
data$ProprieteSol.Sable_.[is.na(data$ProprieteSol.Sable_.)] = data$SableCentroid[is.na(data$ProprieteSol.Sable_.)]
data$ProprieteSol.Limon_. = 100 - data$ProprieteSol.Argile_. - data$ProprieteSol.Sable_.


########
# CoDa #
########
# Texture + C
soilTextParts = cbind(data$ProprieteSol.C_._imp,
                      (100 - data$ProprieteSol.C_._imp) * data[c('ProprieteSol.Sable_.',
                                                                 'ProprieteSol.Limon_.',
                                                                 'ProprieteSol.Argile_.')] / 100)
rowSums(soilTextParts)

soilTextComp = acomp(soilTextParts)
soilTextSBP = matrix(c(-1, 1, 1, 1,
                       0, 1, 1,-1,
                       0, 1,-1, 0),
                     byrow = TRUE,
                     ncol = 4)
colnames(soilTextSBP) = c('C', 'S', 's', 'c') # C: Carbone (carbon), S: Sable (sand), s: Limon (silt), c: Argile (clay)
soilTextBal = ilr(soilTextComp, V = gsi.buildilrBase(t(soilTextSBP)))
colnames(soilTextBal) = c("[C| c,s,S]", "[c | s,S]", "[s | S]")
soilTextBal = ilrNA(soilTextComp, soilTextSBP, soilTextBal)

data$Texture_C.asS = soilTextBal[, 1]
data$Texture_a.sS = soilTextBal[, 2]
data$Texture_s.S = soilTextBal[, 3]


##########################################
# Mise à jour des variables explicatives #
##########################################
variable.explicative

# Changement de nom
variable.explicative[variable.explicative == "ProprieteSol.pHeauCaCl2"] = "pH"
variable.explicative[variable.explicative == "ProprieteSol.C_."] = "Texture_C.asS"
variable.explicative[variable.explicative == "ProprieteSol.Argile_."] = "Texture_a.sS"
variable.explicative[variable.explicative == "ProprieteSol.Sable_."] = "Texture_s.S"

# Ajouts
variable.explicative.bal = c(variable.explicative,
                         'Sol_Al.P', 'Sol_Fv.AlP',
                         'SerieSol_G.P', 'SerieSol_LG.SG')


#################
# Reshape table #
#################
library(tidyr)
data_rs = tbl_df(data[c('ID',
                        variable.explicative.bal,
                        selection.climat,
                        selection.climat_.30,
                        variable.performance,
                        variables.dose)])
#data_rs %>% print(n = Inf)
data.frame(names(data_rs))

data_rs = data_rs %>%
  gather(traitementP, doseP, which(grepl(x = colnames(data_rs), pattern = "Dose")))
data_rs = data_rs[order(data_rs$ID, data_rs$doseP), ]
data_rs$traitementP
data.frame(names(data_rs))

variable.performance.generic = c(
  "Performance.P.Rendement_t.ha",
  "Performance.P.DensiteGrain_g.L",
  "Performance.P.HumiditeGrain_.",
  "Performance.P.HauteurPlant_cm"
)

for (i in seq_along(variable.performance.generic)) {
  ndose = data_rs %>%
    select(starts_with(variable.performance.generic[i])) %>%
    ncol()
  types = c('temoin', paste0('dose', 1:(ndose-1)))
  data_rs[variable.performance.generic[i]] = NA
  for (j in types) {
    data_rs[data_rs$traitementP == paste0('Dose.P.Quantite_kgP2o5.ha_', j), variable.performance.generic[i]] =
      data_rs[data_rs$traitementP == paste0('Dose.P.Quantite_kgP2o5.ha_', j),
              paste0(variable.performance.generic[i], '_', j)]
  }
}

data_rs = data_rs[c('ID',
                    variable.explicative.bal,
                    selection.climat,
                    selection.climat_.30,
                    variable.performance.generic,
                    'traitementP',
                    'doseP')]
data_rs$ID
data.frame(data_rs$doseP, data_rs$Performance.P.Rendement_t.ha)


###################################################################################
# Vétifier l'influence de hauteur du plant au démarrage sur la performance finale #
###################################################################################

# Calcul des performances relatives (rendement, densité des grains, humiddit des grains)
# sous forme de log ratio relativement au témoin, log10(performance.dose / performance.témoin)
data_rs$rdt_log.temoin = rep(NA, nrow(data_rs))
data_rs$dg_log.temoin = rep(NA, nrow(data_rs))
data_rs$hg_log.temoin = rep(NA, nrow(data_rs))
data_rs$hp_log.temoin = rep(NA, nrow(data_rs))
id.fac = factor(data_rs$ID)
for (i in 1:nlevels(id.fac)) {
  # rendement
  rdt.ref = data$Performance.P.Rendement_t.ha_temoin[data$ID == levels(id.fac)[i]]
  rdt.dose = data_rs$Performance.P.Rendement_t.ha[data_rs$ID == levels(id.fac)[i]]
  data_rs$rdt_log.temoin[data_rs$ID == levels(id.fac)[i]] = log10(rdt.dose / rdt.ref)
  # densité des grains
  dg.ref = data$Performance.P.DensiteGrain_g.L_temoin[data$ID == levels(id.fac)[i]]
  dg.dose = data_rs$Performance.P.DensiteGrain_g.L[data_rs$ID == levels(id.fac)[i]]
  data_rs$dg_log.temoin[data_rs$ID == levels(id.fac)[i]] = log10(dg.dose / dg.ref)
  # humidité des grains
  hg.ref = data$Performance.P.HumiditeGrain_._temoin[data$ID == levels(id.fac)[i]]
  hg.dose = data_rs$Performance.P.HumiditeGrain_.[data_rs$ID == levels(id.fac)[i]]
  data_rs$hg_log.temoin[data_rs$ID == levels(id.fac)[i]] = log10(hg.dose / hg.ref)
  # hauteur des plans
  hp.ref = data$Performance.P.HauteurPlant_cm_temoin[data$ID == levels(id.fac)[i]]
  hp.dose = data_rs$Performance.P.HauteurPlant_cm[data_rs$ID == levels(id.fac)[i]]
  data_rs$hp_log.temoin[data_rs$ID == levels(id.fac)[i]] = log10(hp.dose / hp.ref)
}

# Modèle linéaire pour évaluer l'influence de la hauteur au démarrage sur le rendement
## Relatif
### Sélection des données
data_mm_preprocess = data_rs[c('ID', 
                               "rdt_log.temoin",
                               'hp_log.temoin')]
data_mm_preprocess$ID = factor(data_mm_preprocess$ID)
data_mm = data_mm_preprocess %>%
  na.omit() %>%
  droplevels()

### Exploration
plot(x = data_mm$hp_log.temoin, y = data_mm$rdt_log.temoin) # exploratoire

### Modèle relatif de rendement
mm_rdtrel_lin.hp = lme(fixed = rdt_log.temoin ~ hp_log.temoin,
                    data = data_mm, 
                    random =  ~ 1 | ID,
                    method = 'REML')
rsq(y = data_mm$rdt_log.temoin,
    y_hat = predict(mm_rdtrel_lin.hp, level = 0))
intervals(mm_rdtrel_lin.hp)
summary(mm_rdtrel_lin.hp)$tTable

## Absolu
### Sélection des données
data_mm_preprocess = data_rs[c('ID', 
                               "Performance.P.Rendement_t.ha",
                               'Performance.P.HauteurPlant_cm')]
data_mm_preprocess$ID = factor(data_mm_preprocess$ID)
data_mm = data_mm_preprocess %>%
  na.omit() %>%
  droplevels()

### Exploration
plot(x = data_mm$Performance.P.HauteurPlant_cm, y = data_mm$Performance.P.Rendement_t.ha,
     xlab = "Hauteur du plant, stade 5-6 feuilles (cm)",
     ylab = "Rendement (t/ha)") 
abline(v = 80, lty = 2)
abline(h = 8, lty = 2)

### Modèle absolu de rendement
mm_rdtabs_lin.hp = lme(fixed = Performance.P.Rendement_t.ha ~ Performance.P.HauteurPlant_cm,
                       data = data_mm, 
                       random =  ~ 1 | ID,
                       method = 'REML')
rsq(y = data_mm$Performance.P.Rendement_t.ha,
    y_hat = predict(mm_rdtabs_lin.hp, level = 0))
intervals(mm_rdtabs_lin.hp)
summary(mm_rdtabs_lin.hp)$tTable

### Résultats graphique
plot(x = data_mm$Performance.P.HauteurPlant_cm, y = data_mm$Performance.P.Rendement_t.ha,
     xlab = "Hauteur du plant, stade 5-6 feuilles (cm)",
     ylab = "Rendement (t/ha)") 
#abline(v = 80, lty = 2)
#abline(h = 8, lty = 2)
#lines(x = data_mm$Performance.P.HauteurPlant_cm,
 #     y = predict(mm_rdtabs_lin.hp, level = 0))

data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.Rendement_t.ha >= 8
sum(data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.Rendement_t.ha >= 8)
nTN = sum(data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.Rendement_t.ha >= 8)
sum(data_mm$Performance.P.HauteurPlant_cm >= 80)
nTN / sum(data_mm$Performance.P.HauteurPlant_cm >= 80)


data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.Rendement_t.ha < 8
sum(data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.Rendement_t.ha < 8)
nTP = sum(data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.Rendement_t.ha < 8)
sum(data_mm$Performance.P.HauteurPlant_cm < 80)
nTP / sum(data_mm$Performance.P.HauteurPlant_cm < 80)



# Modèle linéaire pour évaluer l'influence de la hauteur au démarrage sur la densité des grains
## Relatif
### Sélection des données
data_mm_preprocess = data_rs[c('ID', 
                               "dg_log.temoin",
                               'hp_log.temoin')]
data_mm_preprocess$ID = factor(data_mm_preprocess$ID)
data_mm = data_mm_preprocess %>%
  na.omit() %>%
  droplevels()

### Exploration
plot(x = data_mm$hp_log.temoin, y = data_mm$dg_log.temoin) # exploratoire

### Modèle relatif densit?
mm_dgrel_lin.hp = lme(fixed = dg_log.temoin ~ hp_log.temoin,
                       data = data_mm, 
                       random =  ~ 1 | ID,
                       method = 'REML')
rsq(y = data_mm$dg_log.temoin,
    y_hat = predict(mm_dgrel_lin.hp, level = 0))
intervals(mm_dgrel_lin.hp)
summary(mm_dgrel_lin.hp)$tTable

## Absolu
### Sélection des données
data_mm_preprocess = data_rs[c('ID', 
                               "Performance.P.DensiteGrain_g.L",
                               'Performance.P.HauteurPlant_cm')]
data_mm_preprocess$ID = factor(data_mm_preprocess$ID)
data_mm = data_mm_preprocess %>%
  na.omit() %>%
  droplevels()

### Exploration
plot(x = data_mm$Performance.P.HauteurPlant_cm, y = data_mm$Performance.P.DensiteGrain_g.L,
     xlab = "Hauteur du plant, stade 5-6 feuilles (cm)",
     ylab = "Densit? du grain (g/L)") 
abline(v = 80, lty = 2)
abline(h = 700, lty = 2)

### Modèle absolu densit?
mm_dgabs_lin.hp = lme(fixed = Performance.P.DensiteGrain_g.L ~ Performance.P.HauteurPlant_cm,
                       data = data_mm, 
                       random =  ~ 1 | ID,
                       method = 'REML')
rsq(y = data_mm$Performance.P.DensiteGrain_g.L,
    y_hat = predict(mm_dgabs_lin.hp, level = 0))
intervals(mm_dgabs_lin.hp)
summary(mm_dgabs_lin.hp)$tTable

### Résultats graphique
plot(x = data_mm$Performance.P.HauteurPlant_cm, y = data_mm$Performance.P.DensiteGrain_g.L,
     xlab = "Hauteur du plant, stade 5-6 feuilles (cm)",
     ylab = "Densit? du grain (g/L)") 
#abline(v = 80, lty = 2)
#abline(h = 675, lty = 2)
lines(x = data_mm$Performance.P.HauteurPlant_cm,
      y = predict(mm_dgabs_lin.hp, level = 0))

data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.DensiteGrain_g.L >= 720
sum(data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.DensiteGrain_g.L >= 720)
nTN = sum(data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.DensiteGrain_g.L >= 720)
sum(data_mm$Performance.P.HauteurPlant_cm >= 80)
nTN / sum(data_mm$Performance.P.HauteurPlant_cm >= 80)
 
data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.DensiteGrain_g.L < 720
sum(data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.DensiteGrain_g.L < 720)
nTP = sum(data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.DensiteGrain_g.L < 720)
sum(data_mm$Performance.P.HauteurPlant_cm < 80)
nTP / sum(data_mm$Performance.P.HauteurPlant_cm < 80)



density_class = cut(data_rs$Performance.P.DensiteGrain_g.L, c(0, 700, 1000), labels = c("low density", "high density"))

for (i in names(data_rs)) {
  if(class(data_rs[[i]]) == "numeric") {
    png(filename = paste0('image_density-class/trial_', i, '.png'),
        width = 500, height = 500, res = 90)
    boxplot(data_rs[[i]] ~ density_class, main = paste0('trial_', i))
    dev.off()
  }
}



# Modèle linéaire pour évaluer l'influence de la hauteur au démarrage sur la humidité des grains
## Relatif
### Sélection des données
data_mm_preprocess = data_rs[c('ID', 
                               "hg_log.temoin",
                               'hp_log.temoin')]
data_mm_preprocess$ID = factor(data_mm_preprocess$ID)
data_mm = data_mm_preprocess %>%
  na.omit() %>%
  droplevels()

### Exploration
plot(x = data_mm$hp_log.temoin, y = data_mm$hg_log.temoin) # exploratoire

### Modèle relatif de rendement
mm_hgrel_lin.hp = lme(fixed = hg_log.temoin ~ hp_log.temoin,
                      data = data_mm, 
                      random =  ~ 1 | ID,
                      method = 'REML')
rsq(y = data_mm$hg_log.temoin,
    y_hat = predict(mm_hgrel_lin.hp, level = 0))
intervals(mm_hgrel_lin.hp)
summary(mm_hgrel_lin.hp)$tTable

## Absolu
### Sélection des données
data_mm_preprocess = data_rs[c('ID', 
                               "Performance.P.HumiditeGrain_.",
                               'Performance.P.HauteurPlant_cm')]
data_mm_preprocess$ID = factor(data_mm_preprocess$ID)
data_mm = data_mm_preprocess %>%
  na.omit() %>%
  droplevels()

### Exploration
plot(x = data_mm$Performance.P.HauteurPlant_cm, y = data_mm$Performance.P.HumiditeGrain_.,
     xlab = "Hauteur du plant, stade 5-6 feuilles (cm)",
     ylab = "Humidit? du grain (%)") 
abline(v = 80, lty = 2)
abline(h = 32, lty = 2)

### Modèle absolu humidit?
mm_hgabs_lin.hp = lme(fixed = Performance.P.HumiditeGrain_. ~ Performance.P.HauteurPlant_cm,
                      data = data_mm, 
                      random =  ~ 1 | ID,
                      method = 'REML')
rsq(y = data_mm$Performance.P.HumiditeGrain_.,
    y_hat = predict(mm_hgabs_lin.hp, level = 0))
intervals(mm_hgabs_lin.hp)
summary(mm_hgabs_lin.hp)$tTable

### Résultats graphique
plot(x = data_mm$Performance.P.HauteurPlant_cm, y = data_mm$Performance.P.HumiditeGrain_.,
     xlab = "Hauteur du plant, stade 5-6 feuilles (cm)",
     ylab = "Humidit? du grain (%)") 
#abline(v = 80, lty = 2)
#abline(h = 32, lty = 2)
lines(x = data_mm$Performance.P.HauteurPlant_cm,
      y = predict(mm_hgabs_lin.hp, level = 0))


data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.HumiditeGrain_. <= 32
sum(data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.HumiditeGrain_. <= 32)
nTN = sum(data_mm$Performance.P.HauteurPlant_cm >= 80 & data_mm$Performance.P.HumiditeGrain_. <= 32)
sum(data_mm$Performance.P.HauteurPlant_cm >= 80)
nTN / sum(data_mm$Performance.P.HauteurPlant_cm >= 80)


data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.HumiditeGrain_. > 32
sum(data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.HumiditeGrain_. > 32)
nTP = sum(data_mm$Performance.P.HauteurPlant_cm < 80 & data_mm$Performance.P.HumiditeGrain_. > 32)
sum(data_mm$Performance.P.HauteurPlant_cm < 80)
nTP / sum(data_mm$Performance.P.HauteurPlant_cm < 80)


##########
# Modèle #
##########
ggplot(filter(data_rs, ID %in% 1:50), aes(x = doseP, y = Performance.P.Rendement_t.ha)) +
  geom_line(aes(colour = factor(ID)))


# Rendement
# ---------

## Préparer les variables explicatives
## Définir var_mm_rdt pour inclure ou exclure des variables
var_mm_rdt = c(variable.explicative.bal[variable.explicative.bal != 'UTM.FA' & 
                                          variable.explicative.bal != 'UTM.diff'],
               selection.climat)
data_mm_rdt_preprocess = data_rs[c('ID', variable.performance.generic, 'doseP', var_mm_rdt)]
par(mar = c(5, 15, 1, 1)) 
barplot(apply(data_mm_rdt_preprocess, 2, function(x) 100 * sum(is.na(x)) / length(x)),
        horiz = TRUE, las = 1, xlab = "Pourcentage de données manquantes")
dev.off()
data_mm_rdt_preprocess$ID = factor(data_mm_rdt_preprocess$ID)

## Filtrer les données pour enlever les valeurs manquantes
data_mm_rdt = data_mm_rdt_preprocess %>%
  select(-Performance.P.DensiteGrain_g.L, -Performance.P.HumiditeGrain_., -Performance.P.HauteurPlant_cm) %>%
  na.omit() %>%
  droplevels()

prep_nlme_rdt = prepare_nlme(X = data_mm_rdt, variable = var_mm_rdt) # [var_mm != 'Culture.Densite.Plan_plant.ha'] # [var_mm != 'SDI']
start_vector_rdt = prep_nlme_rdt$start
rhs_rdt = prep_nlme_rdt$rhs
data_mm_rdt = prep_nlme_rdt$scaled

## Modèle Mitcherlcih
mm_rdt = nlme(Performance.P.Rendement_t.ha ~ Asym * ((1 - exp(-Rate * (doseP + Envi)))),
              data = data_mm_rdt, 
              start = c(Asym = c(9, start_vector_rdt),
                        Rate = c(0.06, start_vector_rdt),
                        Envi = c(80, start_vector_rdt) 
              ), 
              fixed = list(as.formula(paste('Asym ~', rhs_rdt)),
                           as.formula(paste('Rate ~', rhs_rdt)),
                           as.formula(paste('Envi ~', rhs_rdt))
              ), 
              random = Asym ~ 1 | ID,
              control = list(maxIter = 100, returnObject=TRUE, 
                             msVerbose=TRUE, minScale = 1e-10),
              method = 'REML')
AIC(mm_rdt)

rsq(y = data_mm_rdt$Performance.P.Rendement_t.ha,
    y_hat = predict(mm_rdt))
rsq(y = data_mm_rdt$Performance.P.Rendement_t.ha,
    y_hat = predict(mm_rdt, level = 0))

var_mm_rdt
dfcat = data.frame(
  varCat = c('.(Intercept)', var_mm_rdt),
  varCatNamesTo = c('Intercept', 'Management', 'Soil chemistry',
                    'Soil texture', 'Soil texture', 'Soil texture',
                    'Soil chemistry', 'Soil chemistry',
                   'Soil pedology', 'Soil pedology', 
                   'Climate', 'Climate', 'Climate', 'Climate'),
  stringsAsFactors = FALSE)

p_mm_rdt = nlmePlotConfTable(mm = mm_rdt,
                             conf_level = 0.95,
                             pval.breaks = c(0, 0.05, 0.1, 1.00), # p-value breaks
                             varCat = dfcat$varCat,
                             varCatNamesTo = dfcat$varCatNamesTo,
                             varCatOrder = 1:length(unique(dfcat$varCatNamesTo)),
                             paramCat = c('Asym', 'Rate', 'Envi'), # name of parameters in the model
                             greyEnd = 0.65,
                             removeIntercept = TRUE,
                             showRanef = TRUE,
                             numApart = FALSE,
                             varImpType = "cat",
                             limits_df = NULL,
                             #modify.param.names = c("Intercept", "Slope"),
                             varNamesFrom = translate_col$from_name_mm,
                             varNamesTo = translate_col$to_name_mm_long_en)
p_mm_rdt[[1]]
p_mm_rdt[[2]]

summary(mm_rdt)$tTable
plot(density(ranef(mm_rdt)[[1]]))
mm_rdt_ranef_ci = t.test(ranef(mm_rdt)[[1]])$conf.int
abline(v = mm_rdt_ranef_ci[1], col = "red")
abline(v = mm_rdt_ranef_ci[2], col = "red")

### Graphiques de Mitshcerlich appliqués sur les essais
#source('fonctions.R')

#### R² at level 1 for all tests
rsqMM <- c()
for (i in 1:length(levels(data_mm_rdt$ID))) {
  ess <- levels(data_mm_rdt$ID)[i]
  rsqMM[i] <- rsq(y = data_mm_rdt$Performance.P.Rendement_t.ha[data_mm_rdt$ID == ess], 
                  y_hat = predict(mm_rdt, level = 1)[data_mm_rdt$ID == ess])
}
data.frame(levels(data_mm_rdt$ID), rsqMM)

trials <- levels(data_mm_rdt$ID)
#trials <- c('4', '71', '89', '342', '400')
m_pred_id = list()

for (i in 1:length(trials)) {
  set.seed(555)
  ess <- trials[i]
  newdata = data_mm_rdt %>% 
    filter(ID == ess) %>%
    select(one_of(var_mm_rdt)) %>%
    summarise_each(funs(median))
  newdata$ID = ess
  casBase = newdata
  nbNd = 50
  for (i in 2:nbNd) {
    newdata <- rbind(newdata, casBase)
  }
  newdata$doseP <- seq(0, 225, length = nbNd)
  m_pred_id[[ess]] = mitschPred(mm = mm_rdt, newdata = newdata, rhs = rhs_rdt, col_dose = 'doseP', 
                         ranEf = ranef(mm_rdt)[rownames(ranef(mm_rdt)) == ess, 1])
  png(filename = paste0('image_mitsch_rdt/trial_', ess, '.png'),
      width = 500, height = 500, res = 90)
  plot(newdata$doseP, m_pred_id[[ess]]$pred, 
       xlim = c(0, max(newdata$doseP)), ylim = c(4, 14),
       main = paste('Trial', ess),
       xlab = 'Dose P (kg/ha)',
       ylab = 'Yield (t/ha)', type = 'l')
  points(x = data_mm_rdt$doseP[data_mm_rdt$ID == ess],
         y = data_mm_rdt$Performance.P.Rendement_t.ha[data_mm_rdt$ID == ess])
  dev.off()
}

## Modèle linéaire
data_mm_rdt$dose.P_sc = (data_mm_rdt$doseP - mean(data_mm_rdt$doseP)) / sd(data_mm_rdt$doseP)
mean(data_mm_rdt$dose.P_sc)
sd(data_mm_rdt$dose.P_sc)

mm_rdt_lin = lme(fixed = as.formula(paste('Performance.P.Rendement_t.ha ~', rhs_rdt,
                                          '+ dose.P_sc')),
                 data = data_mm_rdt, 
                 random =  ~ 1 | ID,
                 method = 'REML')

rsq(y = data_mm_rdt$Performance.P.Rendement_t.ha,
    y_hat = predict(mm_rdt_lin))
rsq(y = data_mm_rdt$Performance.P.Rendement_t.ha,
    y_hat = predict(mm_rdt_lin, level = 0))

dfcat = data.frame(
  varCat = c('.(Intercept)', var_mm_rdt),
  varCatNamesTo = c('Intercept', 'Management', 'Soil chemistry',
                    'Soil texture', 'Soil texture', 'Soil texture',
                    'Soil chemistry', 'Soil chemistry',
                    'Soil pedology', 'Soil pedology', 
                    'Climate', 'Climate', 'Climate', 'Climate'),
  stringsAsFactors = FALSE)

p_lmm_rdt = nlmePlotConfTable(mm = mm_rdt_lin,
                              conf_level = 0.95,
                              pval.breaks = c(0, 0.05, 0.1, 1.00), # p-value breaks
                              varCat = dfcat$varCat,
                              varCatNamesTo = dfcat$varCatNamesTo,
                              varCatOrder = 1:length(unique(dfcat$varCatNamesTo)),
                              paramCat = c(NA), # name of parameters in the model
                              greyEnd = 0.65,
                              removeIntercept = TRUE,
                              showRanef = TRUE,
                              numApart = FALSE,
                              varImpType = "cat",
                              limits_df = NULL,
                              modify.param.names = NULL, #c("Intercept", "Slope"),
                              varNamesFrom = translate_col$from_name_mm,
                              varNamesTo = translate_col$to_name_mm_long_en)
p_lmm_rdt[[1]]
p_lmm_rdt[[2]]
AIC(mm_rdt_lin)

## Exporter les modèles en image

trials <- levels(data_mm_rdt$ID)
m_pred_id = list()
for (i in 1:length(trials)) {
  ess <- trials[i]
  newdata = data_mm_rdt %>% 
    filter(ID == ess) %>%
    select(one_of(var_mm_rdt)) %>%
    summarise_each(funs(median))
  newdata$ID = ess
  casBase = newdata
  nbNd = 3
  for (i in 2:nbNd) {
    newdata <- rbind(newdata, casBase)
  }
  newdata$dose.P_sc <- seq(-1.1, 2, length = nbNd)
  newdata$doseP = newdata$dose.P_sc * sd(data_mm_rdt$doseP) + mean(data_mm_rdt$doseP) # échelle originale
  m_pred_id[[ess]] = predict(mm_rdt, newdata = newdata)
  png(filename = paste0('image_lin_rdt/trial_', ess, '.png'),
      width = 500, height = 500, res = 90)
  plot(newdata$doseP, m_pred_id[[ess]], 
       xlim = c(0, max(newdata$doseP)),
       ylim = c(4, 16),
       main = paste('Trial', ess),
       xlab = 'Dose P (kg/ha)',
       ylab = 'RENDEMENT (t/ha)', type = 'l')
  points(x = data_mm_rdt$doseP[data_mm_rdt$ID == ess],
         y = data_mm_rdt$Performance.P.Rendement_t.ha[data_mm_rdt$ID == ess])
  dev.off()
}



##### NOTE rdt #####
# Sélectionner Mitscherlich
# car il fournit beaucoup d'information sur l'influence 
# des variables sur les 3 paramètres
# Présenter aussi le modèle linéaire
# On y voit que la dose a peu (mais de manière significative) d'influence
# su le rendement. Contrairement aux balances [Fv | Al,P] et [Al | P], qui eux ont 
# une influence majeure et significative. Ainsi, la dose de P doit être analysée comme
# élément permettant de soutenir indirectement le rendement en soutenant les balances
# [Fv | Al,P] et [Al | P]


# Densité des grains
# -----------------

## Préparer les variables explicatives
## Définir var_mm_rdt pour inclure ou exclure des variables
var_mm_dg = c(variable.explicative.bal[variable.explicative.bal != 'UTM.FA' & 
                                          variable.explicative.bal != 'UTM.diff'],
               selection.climat)
data_mm_dg_preprocess = data_rs[c('ID', variable.performance.generic, 'doseP', var_mm_dg)]
data_mm_dg_preprocess$ID = factor(data_mm_dg_preprocess$ID)

## Filtrer les données pour enlever les valeurs manquantes
data_mm_dg = data_mm_dg_preprocess %>%
  select(-Performance.P.Rendement_t.ha, -Performance.P.HumiditeGrain_., -Performance.P.HauteurPlant_cm) %>%
  na.omit() %>%
  droplevels()

prep_nlme_dg = prepare_nlme(X = data_mm_dg, variable = var_mm_dg)
start_vector_dg = prep_nlme_dg$start
rhs_dg = prep_nlme_dg$rhs
data_mm_dg = prep_nlme_dg$scaled

## Modèle Mitscherlich
mm_dg = nlme(Performance.P.DensiteGrain_g.L ~ Asym * ((1 - exp(-Rate * (doseP + Envi)))),
             data = data_mm_dg, 
             start = c(Asym = c(700, start_vector_dg),
                       Rate = c(0.05, start_vector_dg),
                       Envi = c(50, start_vector_dg) 
             ), 
             fixed = list(as.formula(paste('Asym ~', rhs_dg)),
                          as.formula(paste('Rate ~', rhs_dg)),
                          as.formula(paste('Envi ~', rhs_dg))
             ), 
             random = Asym ~ 1 | ID,
             control = list(maxIter = 50, returnObject=TRUE, 
                            msVerbose=TRUE, minScale = 1e-8),
             method = 'REML')

rsq(y = data_mm_dg$Performance.P.DensiteGrain_g.L,
    y_hat = predict(mm_dg))
rsq(y = data_mm_dg$Performance.P.DensiteGrain_g.L,
    y_hat = predict(mm_dg, level = 0))


dfcat = data.frame(
  varCat = c('.(Intercept)', var_mm_dg),
  varCatNamesTo = c('Intercept', 'Management', 'Soil chemistry',
                    'Soil texture', 'Soil texture', 'Soil texture',
                    'Soil chemistry', 'Soil chemistry',
                    'Soil pedology', 'Soil pedology', 
                    'Climate', 'Climate', 'Climate', 'Climate'),
  stringsAsFactors = FALSE)

p_mm_dg = nlmePlotConfTable(mm = mm_dg,
                             conf_level = 0.95,
                             pval.breaks = c(0, 0.05, 0.1, 1.00),
                             varCat = dfcat$varCat,
                             varCatNamesTo = dfcat$varCatNamesTo,
                             varCatOrder = 1:length(unique(dfcat$varCatNamesTo)),
                            paramCat = c('Asym', 'Rate', 'Envi'), # name of parameters in the model
                            greyEnd = 0.65,
                            removeIntercept = TRUE,
                            showRanef = TRUE,
                            numApart = FALSE,
                            varImpType = "cat",
                            limits_df = NULL,
                            #modify.param.names = c("Intercept", "Slope"),
                            varNamesFrom = translate_col$from_name_mm,
                            varNamesTo = translate_col$to_name_mm_long_en)
p_mm_dg[[1]]
p_mm_dg[[2]]

## Modèle linéaire
data_mm_dg$dose.P_sc = (data_mm_dg$doseP - mean(data_mm_dg$doseP)) / sd(data_mm_dg$doseP)

mm_dg_lin = lme(fixed = as.formula(paste('Performance.P.DensiteGrain_g.L ~', rhs_dg, '+ dose.P_sc')),
                data = data_mm_dg, 
                random =  ~ 1 | ID,
                method = 'REML')

summary(mm_dg_lin)$tTable

rsq(y = data_mm_dg$Performance.P.DensiteGrain_g.L,
    y_hat = predict(mm_dg_lin))
rsq(y = data_mm_dg$Performance.P.DensiteGrain_g.L,
    y_hat = predict(mm_dg_lin, level = 0))
p_lmm_dg = nlmePlotConfTable(mm = mm_dg_lin,
                            conf_level = 0.95,
                            pval.breaks = c(0, 0.05, 0.1, 1.00),
                            varCat = dfcat$varCat,
                            varCatNamesTo = dfcat$varCatNamesTo,
                            varCatOrder = 1:length(unique(dfcat$varCatNamesTo)),
                             paramCat = c(NA), # name of parameters in the model
                             greyEnd = 0.65,
                             removeIntercept = TRUE,
                             showRanef = TRUE,
                             numApart = FALSE,
                             varImpType = "cat",
                             limits_df = NULL,
                             modify.param.names = NULL, #c("Intercept", "Slope"),
                             varNamesFrom = translate_col$from_name_mm,
                             varNamesTo = translate_col$to_name_mm_long_en)
p_lmm_dg[[1]]
p_lmm_dg[[2]]

## Exporter les modèles en image
trials <- levels(data_mm_dg$ID)
m_pred_id = list()
for (i in 1:length(trials)) {
  ess <- trials[i]
  newdata = data_mm_dg %>% 
    filter(ID == ess) %>%
    select(one_of(var_mm_dg)) %>%
    summarise_each(funs(median))
  newdata$ID = ess
  casBase = newdata
  nbNd = 3
  for (i in 2:nbNd) {
    newdata <- rbind(newdata, casBase)
  }
  newdata$dose.P_sc <- seq(-1.1, 2, length = nbNd)
  newdata$doseP = newdata$dose.P_sc * sd(data_mm_dg$doseP) + mean(data_mm_dg$doseP) # échelle originale
  m_pred_id[[ess]] = predict(mm_dg_lin, newdata = newdata)
  png(filename = paste0('image_lin_dg/trial_', ess, '.png'),
      width = 500, height = 500, res = 90)
  plot(newdata$doseP, m_pred_id[[ess]], 
       xlim = c(0, max(newdata$doseP)),
       ylim = c(520, 825),
       main = paste('Trial', ess),
       xlab = 'Dose P (kg/ha)',
       ylab = 'DENSITE DU GRAIN (g/l)', type = 'l')
  points(x = data_mm_dg$doseP[data_mm_dg$ID == ess],
         y = data_mm_dg$Performance.P.DensiteGrain_g.L[data_mm_dg$ID == ess])
  dev.off()
}


##### NOTE dg #####
# Sélectionnerle modèle linéaire
# car le modèle de Mitscherlich ne converge pas, et peu d'information semble s'en dégager
# le modèle linéaire permet de mieux apprécier l'influence des variables sur la
# densité des grains, et de comparer leur influence à celui de la dose de P.


# Humidité des grains
#--------------------
var_mm_hg = c(variable.explicative.bal[variable.explicative.bal != 'UTM.FA' & 
                                         variable.explicative.bal != 'UTM.diff'],
              selection.climat)
data_mm_hg_preprocess = data_rs[c('ID', variable.performance.generic, 'doseP', var_mm_hg)]
data_mm_hg_preprocess$ID = factor(data_mm_hg_preprocess$ID)

## Filtrer les données pour enlever les valeurs manquantes
data_mm_hg = data_mm_hg_preprocess %>%
  select(-Performance.P.Rendement_t.ha, -Performance.P.Rendement_t.ha, -Performance.P.HauteurPlant_cm) %>%
  na.omit() %>%
  droplevels()

prep_nlme_hg = prepare_nlme(X = data_mm_hg, variable = var_mm_hg)
start_vector_hg = prep_nlme_hg$start
rhs_hg = prep_nlme_hg$rhs
data_mm_hg = prep_nlme_hg$scaled

## Modèle Mitscherlich
mm_hg = nlme(Performance.P.HumiditeGrain_. ~ Asym * ((1 - exp(-Rate * (doseP + Envi)))),
             data = data_mm_hg, 
             start = c(Asym = c(30, start_vector_hg),
                       Rate = c(0.05, start_vector_hg),
                       Envi = c(50, start_vector_hg) 
             ), 
             fixed = list(as.formula(paste('Asym ~', rhs_hg)),
                          as.formula(paste('Rate ~', rhs_hg)),
                          as.formula(paste('Envi ~', rhs_hg))
             ), 
             random = Asym ~ 1 | ID,
             control = list(maxIter = 50, returnObject=TRUE, 
                            msVerbose=TRUE, minScale = 1e-8),
             method = 'REML')

### Résultats
rsq(y = data_mm_hg$Performance.P.HumiditeGrain_.,
    y_hat = predict(mm_hg))
rsq(y = data_mm_hg$Performance.P.HumiditeGrain_.,
    y_hat = predict(mm_hg, level = 0))

dfcat = data.frame(
  varCat = c('.(Intercept)', var_mm_hg),
  varCatNamesTo = c('Intercept', 'Management', 'Soil chemistry',
                    'Soil texture', 'Soil texture', 'Soil texture',
                    'Soil chemistry', 'Soil chemistry',
                    'Soil pedology', 'Soil pedology', 
                    'Climate', 'Climate', 'Climate', 'Climate'),
  stringsAsFactors = FALSE)

p_mm_hg = nlmePlotConfTable(mm = mm_dg,
                            conf_level = 0.95,
                            pval.breaks = c(0, 0.05, 0.1, 1.00),
                            varCat = dfcat$varCat,
                            varCatNamesTo = dfcat$varCatNamesTo,
                            varCatOrder = 1:length(unique(dfcat$varCatNamesTo)),
                            paramCat = c('Asym', 'Rate', 'Envi'), # name of parameters in the model                            greyEnd = 0.65,
                            removeIntercept = TRUE,
                            showRanef = TRUE,
                            numApart = FALSE,
                            varImpType = "cat",
                            limits_df = NULL,
                            varNamesFrom = translate_col$from_name_mm,
                            varNamesTo = translate_col$to_name_mm_long_en)
p_mm_hg[[1]]

## Modèle linéaire
data_mm_hg$dose.P_sc = (data_mm_hg$doseP - mean(data_mm_hg$doseP)) / sd(data_mm_hg$doseP)

mm_hg_lin = lme(fixed = as.formula(paste('Performance.P.HumiditeGrain_. ~', rhs_hg, '+ dose.P_sc')),
                data = data_mm_hg, 
                random =  ~ 1 | ID,
                method = 'REML')
AIC(mm_hg_lin)


summary(mm_hg_lin)$tTable
intervals(mm_hg_lin)

### Résultats
rsq(y = data_mm_hg$Performance.P.HumiditeGrain_.,
    y_hat = predict(mm_hg_lin))
rsq(y = data_mm_hg$Performance.P.HumiditeGrain_.,
    y_hat = predict(mm_hg_lin, level = 0))

dfcat = data.frame(
  varCat = c('.(Intercept)', var_mm_hg),
  varCatNamesTo = c('Intercept', 'Management', 'Soil chemistry',
                    'Soil texture', 'Soil texture', 'Soil texture',
                    'Soil chemistry', 'Soil chemistry',
                    'Soil pedology', 'Soil pedology', 
                    'Climate', 'Climate', 'Climate', 'Climate'),
  stringsAsFactors = FALSE)

p_lmm_hg = nlmePlotConfTable(mm = mm_hg_lin,
                              conf_level = 0.95,
                              pval.breaks = c(0, 0.05, 0.1, 1.00), # p-value breaks
                              varCat = dfcat$varCat,
                              varCatNamesTo = dfcat$varCatNamesTo,
                              varCatOrder = 1:(length(unique(dfcat$varCatNamesTo))+1),
                              paramCat = c(NA), # name of parameters in the model
                              greyEnd = 0.65,
                              removeIntercept = TRUE,
                              showRanef = TRUE,
                              numApart = FALSE,
                              varImpType = "cat",
                              limits_df = NULL,
                              modify.param.names = NULL, #c("Intercept", "Slope"),
                              varNamesFrom = translate_col$from_name_mm,
                              varNamesTo = translate_col$to_name_mm_long_en)
p_lmm_hg[[1]]



## Exporter les modèles en image
trials <- levels(data_mm_hg$ID)
m_pred_id = list()
for (i in 1:length(trials)) {
  ess <- trials[i]
  newdata = data_mm_hg %>% 
    filter(ID == ess) %>%
    select(one_of(var_mm_hg)) %>%
    summarise_each(funs(median))
  newdata$ID = ess
  casBase = newdata
  nbNd = 3
  for (i in 2:nbNd) {
    newdata <- rbind(newdata, casBase)
  }
  newdata$dose.P_sc <- seq(-1.1, 2, length = nbNd)
  newdata$doseP = newdata$dose.P_sc * sd(data_mm_hg$doseP) + mean(data_mm_hg$doseP) # échelle originale
  m_pred_id[[ess]] = predict(mm_hg_lin, newdata = newdata)
  png(filename = paste0('image_lin_hg/trial_', ess, '.png'),
      width = 500, height = 500, res = 90)
  plot(newdata$doseP, m_pred_id[[ess]], 
       xlim = c(0, max(newdata$doseP)),
       ylim = c(8, 50),
       main = paste('Trial', ess),
       xlab = 'Dose P (kg/ha)',
       ylab = 'HUMIDITE DU GRAIN (%)', type = 'l')
  points(x = data_mm_hg$doseP[data_mm_dg$ID == ess],
         y = data_mm_hg$Performance.P.HumiditeGrain_.[data_mm_hg$ID == ess])
  dev.off()
}





# ---------------------- CORRIGÉ JUSQU'ICI ----------------------------------- #


####################
# Mod?le D?marrage #VOIR LES FACTEURS QUI PEUVENT INFLUENCER LA HAUTEUR DU PLANT AU DEMARRAGE
####################
# DonnÃ©es

var_mm_.30 = c(variable.explicative, selection.climat_.30)
var_mm_.30 = var_mm_.30[var_mm_.30 != 'Culture'] # enleer la culture: tous dans la catÃ©gorie MaÃ¯s grain
var_mm_.30 = c(var_mm_.30[var_mm_.30 != 'UTM.FA' & 
                                var_mm_.30 != 'UTM.diff'  & 
                            var_mm_.30 != 'UTMc' ])
data_mm_preprocess_.30 = data_rs[c('ID', variable.performance.generic, 'doseP', var_mm_.30)]
names(data_mm_preprocess_.30)
data_mm_preprocess_miss_.30 = data_mm_preprocess_.30
names(data_mm_preprocess_miss_.30) = c('ID', 'rdt', 'dg', 'hg', 'hp', 'doseP',
                                       'dens.sem', 'Al|P', 'Fv|AlP', 'G|P',
                                       'LG|SG', 'pH', 'C|asS', 'a|sS', 's|S',
                                       'SDI_.30', 'UTM_.30', 'PCP_.30', 'Tmean_.30') # 'Culture',
aggr(data_mm_preprocess_miss_.30)

data_mm_preprocess_.30$ID = factor(data_mm_preprocess_.30$ID)

# Hauteur du plant
## Filtrer les donnÃ©es pour enlever les valeurs manquantes
data_mm = data_mm_preprocess_.30 %>%
  select(-Performance.P.DensiteGrain_g.L, -Performance.P.HumiditeGrain_., -Performance.P.Rendement_t.ha) %>%
  na.omit() %>%
  droplevels()
nrow(data_mm)
nrow(data_mm_preprocess)

prep_nlme = prepare_nlme(X = data_mm, variable = var_mm_.30) # [var_mm != 'Culture.Densite.Plan_plant.ha'] # [var_mm != 'SDI']
start_vector = prep_nlme$start
rhs = prep_nlme$rhs
data_mm = prep_nlme$scaled

## ModÃ¨le Mitcherlcih
quantile(data_mm$Performance.P.HauteurPlant_cm) # pour choisir l'Asymptote
mm_hp = nlme(Performance.P.HauteurPlant_cm ~ Asym * ((1 - exp(-Rate * (doseP + Envi)))),
             data = data_mm, 
             start = c(Asym = c(100, start_vector),
                       Rate = c(0.05, start_vector),
                       Envi = c(80, start_vector) 
             ), 
             fixed = list(as.formula(paste('Asym ~', rhs)),
                          as.formula(paste('Rate ~', rhs)),
                          as.formula(paste('Envi ~', rhs))
             ), 
             random = Asym ~ 1 | ID,
             control = list(maxIter = 100, returnObject = TRUE, 
                            msVerbose = TRUE, minScale = 1e-10),
             method = 'REML')
AIC(mm_hp)
summary(mm_hp)$tTable

rsq(y = data_mm$Performance.P.HauteurPlant_cm,
    y_hat = predict(mm_hp))
rsq(y = data_mm$Performance.P.HauteurPlant_cm,
    y_hat = predict(mm_hp, level = 0))

translate_col = read.csv(file = "donnees/translate_colnames.csv", header = TRUE, sep = ";")
p_mm_hp = nlmePlotConfTable(mm = mm_hp,
                            conf_level = 0.95,
                            pval.breaks = c(0, 0.05, 0.1, 1.00), # p-value breaks
                            varCat = c('.(Intercept)', 'Sol_', 'pH', 'SerieSol', 'Texture', 
                                       'SDI', 'UTMc', 'PcpCum', 'Tmean',
                                       'Culture.'), # variables which are categorical (Intercept)
                            varCatNamesTo = c('Intercept', 'Soil chemistry', 'Soil chemistry', 'Soil pedology', 'Soil texture',
                                              'Climate', 'Climate', 'Climate', 'Climate', 'Management'),
                            varCatOrder = 1:6, # lenght is nlevels of varCatNamesTo
                            paramCat = c('Asym', 'Rate', 'Envi'), # name of parameters in the model
                            greyEnd = 0.65,
                            removeIntercept = TRUE,
                            showRanef = TRUE,
                            numApart = FALSE,
                            varImpType = "cat",
                            limits_df = NULL,
                            #modify.param.names = c("Intercept", "Slope"),
                            varNamesFrom = translate_col$from_name_mm,
                            varNamesTo = translate_col$to_name_mm_long_en)
p_mm_hp[[1]]


## ModÃ¨le linÃ©aire
data_mm$dose.P_sc = (data_mm$doseP - mean(data_mm$doseP)) / sd(data_mm$doseP)
mean(data_mm$dose.P_sc)
sd(data_mm$dose.P_sc)

mm_hp_lin = lme(fixed = as.formula(paste('Performance.P.HauteurPlant_cm ~', rhs, 
                                         '+ dose.P_sc')),
                data = data_mm, 
                random =  ~ 1 | ID,
                method = 'REML')

rsq(y = data_mm$Performance.P.HauteurPlant_cm,
    y_hat = predict(mm_hp_lin))
rsq(y = data_mm$Performance.P.HauteurPlant_cm,
    y_hat = predict(mm_hp_lin, level = 0))

translate_col = read.csv(file = "donnees/translate_colnames.csv", header = TRUE, sep = ";")

p_lmm_hp = nlmePlotConfTable(mm = mm_hp_lin,
                             conf_level = 0.95,
                             pval.breaks = c(0, 0.05, 0.1, 1.00), # p-value breaks
                             varCat = c('.(Intercept)', 'Sol_', 'pH', 'SerieSol', 'Texture', 
                                        'SDI', 'UTMc', 'PcpCum', 'Tmean',
                                        'Culture.', 'dose.P_sc'), # variables which are categorical (Intercept)
                             varCatNamesTo = c('Intercept', 'Soil chemistry', 'Soil chemistry', 'Soil pedology', 'Soil texture',
                                               'Climate', 'Climate', 'Climate', 'Climate', 'Management', 'Management'),
                             varCatOrder = 1:6, # lenght is nlevels of varCatNamesTo
                             paramCat = c(NA), # name of parameters in the model
                             greyEnd = 0.65,
                             removeIntercept = TRUE,
                             showRanef = TRUE,
                             numApart = FALSE,
                             varImpType = "cat",
                             limits_df = NULL,
                             modify.param.names = NULL, #c("Intercept", "Slope"),
                             varNamesFrom = translate_col$from_name_mm,
                             varNamesTo = translate_col$to_name_mm_long_en)
p_lmm_hp[[1]]


summary(mm_hp)$tTable


