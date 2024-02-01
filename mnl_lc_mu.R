#selecionando diretorio
setwd("C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Modelos\\OD\\Banco de dados")

library(apollo)
library(tidyverse)

rm(list = ls())

parallel::detectCores() #cores available

#ler o banco de dados quando tiver o final
database <- read.csv(file = paste0(getwd(),'\\cleaning\\database_final.csv'))

###Preparing environment
apollo_initialise()

###Options controlling the running of the code
apollo_control = list(
  modelName  ="mnl_lc_2mu 13",
  modelDescr = "2 classes, RUM e muRRM, parametros LOS, one mu parameter and sociodemographics in the class allocation model",
  indivID    ="id_pess",
  mixing     = F,
  nCores     = 10
)

#parametros a serem estimados
apollo_beta <- c(asc1 = 0,
                 asc2 = 0,
                 asc3 = 0,
                 #class 1
                 b_ttauto_1 = -0.15,
                 b_ttbus_1 = -0.04,
                 b_ttmetro_1 = -0.023,
                 b_co_1 = -0.60,
                 #class 2
                 b_ttauto_2 = -0.0578,
                 b_ttbus_2 = -0.052,
                 b_ttmetro_2 = -0.036,
                 b_co_2 = -0.16314,
                 b_mu = 0.2,
                 #class membership
                 s_1 = 0,
                 s_2 = 0,
                 b_pico_al = 0,
                 b_idade_al = 0,
                 b_grau1_al = 0,
                 b_grau2_al = 0,
                 b_sex_al = 0,
                 b_inc_al = 0,
                 b_nomorad_al = 0,
                 b_totviag_al = 0,
                 #sociodemographics
                 b_pico = 0,
                 b_idade= 0,
                 b_grau1 = 0,
                 b_grau2 = 0,
                 b_sex_bus = 0,
                 b_sex_metro = 0,
                 b_inc_bus = 0,
                 b_inc_metro = 0,
                 b_nomorad_bus = 0,
                 b_nomorad_metro = 0,
                 b_totviag_bus = 0,
                 b_totviag_metro = 0
                 )

#parametros fixos
apollo_fixed <- c('asc1','s_1',
                  'b_pico',
                  'b_idade',
                  'b_grau1',
                  'b_grau2',
                  'b_sex_bus',
                  'b_sex_metro',
                  'b_inc_bus',
                  'b_inc_metro',
                  'b_nomorad_bus',
                  'b_nomorad_metro',
                  'b_totviag_bus',
                  'b_totviag_metro'
)

#lendo parametros de outro modelo
setwd('C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\JTBS\\modelo\\resultados\\mnl_lc_mu')

apollo_beta <- apollo_readBeta(apollo_beta,apollo_fixed,'mnl_lc_2mu 12',overwriteFixed = FALSE)

###Grouping latent class parameters
apollo_lcPars = function(apollo_beta, apollo_inputs){
  lcpars = list()
  
  lcpars[['b_ttauto']] = list(b_ttauto_1,b_ttauto_2)
  
  lcpars[['b_ttbus']] = list(b_ttbus_1,b_ttbus_2)
  
  lcpars[['b_ttmetro']] = list(b_ttmetro_1,b_ttmetro_2)
  
  lcpars[['b_co']] = list(b_co_1,b_co_2)
  
   
  ###Class membership probabilities based on s_1, s_2: use of MNL fomrula
  V=list()
  
  V[["class_1"]] = s_1
  V[["class_2"]] = s_2 + b_idade_al*idade + b_pico_al*pico_manha +
    b_sex_al*(sexo == 2) + b_inc_al*lnrenda + b_nomorad_al*no_morad + b_totviag_al*tot_viag + b_grau1_al*grau1 + b_grau2_al*grau2
  
  ###Settings for class membership probabilities 
  mnl_settings = list(
    alternatives = c(class_1=1, class_2=2), 
    avail        = 1, 
    choiceVar    = NA, ###No choice variable as only the formula of MNL is used
    V            = V
  )
  ###Class membership probabilities
  lcpars[["pi_values"]] = apollo_mnl(mnl_settings, functionality="raw")
  return(lcpars)
}

###busca pelos dados de entrada necessarios 
apollo_inputs = apollo_validateInputs()

#cosntruindo a funcao de verossimilhanca
apollo_probabilities <-  function(apollo_beta,apollo_inputs,functionality = 'estimate'){
  
  ###Attaches parameters and data so that variables can be referred by name
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  
  ###Create list for probabilities
  P = list()

  ### Define settings for MNL model component
  mnl_settings = list(
    alternatives = c(auto=3, bus=2, metro=1),
    avail        = list(auto=1, bus=1, metro=1),
    choiceVar    = choice
   )
  
  ### Compute class-specific utilities
  V=list()
  
  V[['auto']]  = asc1 + b_ttauto[[1]]*tt_auto + b_co[[1]]*co_auto + b_pico*pico_manha
  
  V[['bus']]  = asc2 + b_ttbus[[1]]*tt_bus_otp + b_co[[1]]*co_bus + b_idade*idade + 
    b_sex_bus*(sexo == 2) + b_inc_bus*lnrenda + b_nomorad_bus*no_morad + b_totviag_bus*tot_viag + b_grau1*grau1 + b_grau2*grau2
  
  V[['metro']] = asc3 + b_ttmetro[[1]]*tt_metro + b_co[[1]]*co_metro  + b_idade*idade + 
    b_sex_metro*(sexo == 2) + b_inc_metro*lnrenda + b_nomorad_metro*no_morad + b_totviag_metro*tot_viag + b_grau1*grau1 + b_grau2*grau2
  
  ###Calculating probabilities based on MNL function for class 1
  mnl_settings$V = V
  P[[1]] = apollo_mnl(mnl_settings, functionality)
  
  ### Compute class-specific regrets
  R=list()
  
  R[['auto']]  =  asc1 + b_mu*(-log(1 + exp((b_ttbus[[2]]*tt_bus_otp - b_ttauto[[2]]*tt_auto)/b_mu)) -
                                    log(1 + exp((b_ttmetro[[2]]*tt_metro - b_ttauto[[2]]*tt_auto)/b_mu)) + b_mu*(-
                                    log(1 + exp((b_co[[2]]*co_metro - b_co[[2]]*co_auto)/b_mu)) - 
                                    log(1 + exp((b_co[[2]]*co_bus - b_co[[2]]*co_auto)/b_mu)))) + b_pico*pico_manha  
  
  
  R[['bus']]  = asc1 + b_mu*(-log(1 + exp((b_ttauto[[2]]*tt_auto - b_ttbus[[2]]*tt_bus_otp)/b_mu)) -
                                   log(1 + exp((b_ttmetro[[2]]*tt_metro - b_ttbus[[2]]*tt_bus_otp)/b_mu)) + b_mu*(-
                                   log(1 + exp((b_co[[2]]*co_auto - b_co[[2]]*co_bus) / b_mu)) - 
                                   log(1 + exp((b_co[[2]]/b_mu)*(co_metro - co_bus))))) + b_idade*idade + b_grau1*grau1 + b_grau2*grau2 + 
    b_sex_bus*(sexo == 2) + b_inc_bus*lnrenda + b_nomorad_bus*no_morad + b_totviag_bus*tot_viag
  
  R[['metro']] = asc1 + b_mu*(- log(1 + exp((b_ttauto[[2]]*tt_auto - b_ttmetro[[2]]*tt_metro)/b_mu)) - 
                                     log(1 + exp((b_ttbus[[2]]*tt_bus_otp - b_ttmetro[[2]]*tt_metro)/b_mu))  + b_mu*(- 
                                     log(1 + exp((b_co[[2]]*co_auto - b_co[[2]]*co_metro)/b_mu)) - 
                                     log(1 + exp((b_co[[2]]/b_mu)*(co_bus - co_metro))))) + b_idade*idade + b_grau1*grau1 + b_grau2*grau2 + 
    b_sex_metro*(sexo == 2) + b_inc_metro*lnrenda + b_nomorad_metro*no_morad + b_totviag_metro*tot_viag 
  
  ###Calculating probabilities based on MNL function for class 2
  mnl_settings$V = R 
  P[[2]] = apollo_mnl(mnl_settings, functionality)
  
  ###Calculating choice probabilities using class membership and conditional probabilities 
  lc_settings   = list(inClassProb = P, classProb=pi_values)
  P[["model"]] = apollo_lc(lc_settings, apollo_inputs, functionality)
  
  ### Prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  
  return(P)
}

### Optional: searching for starting value
apollo_beta = apollo_searchStart(apollo_beta,
                                 apollo_fixed,
                                 apollo_probabilities,
                                 apollo_inputs,
                                 searchStart_settings=list(nCandidates=20))

# ################################################################# #
#### Estimacao do modelo                                        ####
# ################################################################# #

model_lc_murrm <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs,
                             estimate_settings = list(maxIterations=250) )

###Display the output to console with p-values

apollo_modelOutput(model_lc_murrm,modelOutput_settings = list(printPVal=T ,printT1=T,
                                                              printDiagnostics = T))

# ----------------------------------------------------------------- #
#---- FORMATTED OUTPUT (TO FILE, using model name)               ----
# ----------------------------------------------------------------- #
setwd('C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\JTBS\\modelo\\resultados\\mnl_lc_mu')

apollo_saveOutput(model_lc_murrm, saveOutput_settings = list(printPVal=T, printT1=T,
                                                       printDiagnostics = T)) #ver saveoutput list para mais configuracoes

#VTTS 
setwd('C:\\Users\\calde\\Google Drive\\Mestrado 2019_Gabriel Caldeira\\Artigos\\JTBS\\modelo\\resultados\\mnl_lc_mu')

model_lc <- apollo_loadModel('mnl_lc_mu 10')

#RUM
vtts_auto <- model_lc$estimate[['b_ttauto_1']]/model_lc$estimate[['b_co_1']]*60
vtts_bus <- model_lc$estimate[['b_ttbus_1']]/model_lc$estimate[['b_co_1']]*60
vtts_metro <- model_lc$estimate[['b_ttmetro_1']]/model_lc$estimate[['b_co_1']]*60

vtts_auto
vtts_bus
vtts_metro

#RRM
dtt_auto <- (-model_lc$estimate[['b_ttauto_2']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1) +
            ((-model_lc$estimate[['b_ttauto_2']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1))

dco_auto <- (-model_lc$estimate[['b_co_2']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_bus - database$co_auto))+1) +
            ((-model_lc$estimate[['b_co_2']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_metro - database$co_auto))+1))

vtts_auto <-  dtt_auto/dco_auto*60

dtt_bus <- (-model_lc$estimate[['b_ttbus_2']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1) +
           ((-model_lc$estimate[['b_ttbus_2']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1))

dco_bus <- (-model_lc$estimate[['b_co_2']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_auto - database$co_bus))+1) +
           ((-model_lc$estimate[['b_co_2']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_metro - database$co_bus))+1))

vtts_bus <- dtt_bus/dco_bus*60

dtt_metro <- (-model_lc$estimate[['b_ttmetro_2']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1) +
             (-model_lc$estimate[['b_ttmetro_2']]/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro))+1))

dco_metro <- (-model_lc$estimate[['b_co_2']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_auto - database$co_metro))+1) +
             ((-model_lc$estimate[['b_co_2']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_bus - database$co_metro))+1))

vtts_metro <- dtt_metro/dco_metro*60

mean(vtts_auto)
mean(vtts_bus)
mean(vtts_metro)

vtts_lcmurrm <- data.frame(vtts_auto,vtts_bus,vtts_metro)

vtts_lcmurrm <- vtts_lcmurrm %>% 
  rename('vtts_auto_lcmurrm' = 'vtts_auto',
         'vtts_bus_lcmurrm' = 'vtts_bus',
         'vtts_metro_lcmurrm' = 'vtts_metro')

vtts_models <- cbind(vtts_models,vtts_lcmurrm)

vtts_lcmurrm_socio <- data.frame(vtts_auto,vtts_bus,vtts_metro)

vtts_lcmurrm_socio <- vtts_lcmurrm_socio %>% 
  rename('vtts_auto_lcmurrm_socio' = 'vtts_auto',
         'vtts_bus_lcmurrm_socio' = 'vtts_bus',
         'vtts_metro_lcmurrm_socio' = 'vtts_metro')

vtts_models <- cbind(vtts_models,vtts_lcmurrm_socio)

##RRM PAPER
dtt_auto <- (-model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1) +
  ((-model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1))

dco_auto <- (-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_bus - database$co_auto))+1) +
  ((-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_metro - database$co_auto))+1))

vtts_auto <-  dtt_auto/dco_auto*60

dtt_bus <- (-model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1) +
  ((-model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1))

dco_bus <- (-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_auto - database$co_bus))+1) +
  ((-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_metro - database$co_bus))+1))

vtts_bus <- dtt_bus/dco_bus*60

dtt_metro <- (-model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1) +
  ((-model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1))

dco_metro <- (-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_auto - database$co_metro))+1) +
  ((-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']]*(database$co_bus - database$co_metro))+1))

vtts_metro <- dtt_metro/dco_metro*60

mean(vtts_auto)
mean(vtts_bus)
mean(vtts_metro)

vtts_lcmurrm_paper <- data.frame(vtts_auto,vtts_bus,vtts_metro)

vtts_lcmurrm_paper <- vtts_lcmurrm_paper %>% 
  rename('vtts_auto_lcmurrm_paper' = 'vtts_auto',
         'vtts_bus_lcmurrm_paper' = 'vtts_bus',
         'vtts_metro_lcmurrm_paper' = 'vtts_metro')

vtts_models <- cbind(vtts_models,vtts_lcmurrm)

vtts_lcmurrm_socio <- data.frame(vtts_auto,vtts_bus,vtts_metro)

vtts_lcmurrm_socio <- vtts_lcmurrm_socio %>% 
  rename('vtts_auto_lcmurrm_socio' = 'vtts_auto',
         'vtts_bus_lcmurrm_socio' = 'vtts_bus',
         'vtts_metro_lcmurrm_socio' = 'vtts_metro')

vtts_models <- cbind(vtts_models,vtts_lcmurrm_socio)

##### Elasticities #######

### Use the estimated model to make predictions
predictions_base_rum = apollo_prediction(model_lc, apollo_probabilities, apollo_inputs,prediction_settings = list(modelComponent = 1))
predictions_base_rrm = apollo_prediction(model_lc, apollo_probabilities, apollo_inputs,prediction_settings = list(modelComponent = 2))

####RUM
#tt auto
elast_ttauto <- (1-predictions_base_rum[,'auto'])*model_lc$estimate[['b_ttauto_1']]*database$tt_auto
aggfv_ttauto <- sum(elast_ttauto*(database$fe_via*predictions_base_rum[,'auto']/sum(database$fe_via*predictions_base_rum[,'auto'])))
aggfi_ttauto <- sum(elast_ttauto*(database$fe_pess*predictions_base_rum[,'auto']/sum(database$fe_pess*predictions_base_rum[,'auto'])))
agg_ttauto <- sum(elast_ttauto*(predictions_base_rum[,'auto']/sum(predictions_base_rum[,'auto'])))

#tt bus
elast_ttbus <- (1-predictions_base_rum[,'bus'])*model_lc$estimate[['b_ttbus_1']]*database$tt_bus_otp
aggfv_ttbus <- sum(elast_ttbus*(database$fe_via*predictions_base_rum[,'bus']/sum(database$fe_via*predictions_base_rum[,'bus'])))
aggfi_ttbus <- sum(elast_ttbus*(database$fe_pess*predictions_base_rum[,'bus']/sum(database$fe_pess*predictions_base_rum[,'bus'])))
agg_ttbus <- sum(elast_ttbus*(predictions_base_rum[,'bus']/sum(predictions_base_rum[,'bus'])))

#tt metro
elast_ttmetro <- (1-predictions_base_rum[,'metro'])*model_lc$estimate[['b_ttmetro_1']]*database$tt_metro
aggfv_ttmetro <- sum(elast_ttmetro*(database$fe_via*predictions_base_rum[,'metro']/sum(database$fe_via*predictions_base_rum[,'metro'])))
aggfi_ttmetro <- sum(elast_ttmetro*(database$fe_pess*predictions_base_rum[,'metro']/sum(database$fe_pess*predictions_base_rum[,'metro'])))
agg_ttmetro <- sum(elast_ttmetro*(predictions_base_rum[,'metro']/sum(predictions_base_rum[,'metro'])))

#co auto
elast_coauto <- (1-predictions_base_rum[,'auto'])*model_lc$estimate[['b_co_1']]*database$co_auto
aggfv_coauto <- sum(elast_coauto*(database$fe_via*predictions_base_rum[,'auto']/sum(database$fe_via*predictions_base_rum[,'auto'])))
aggfi_coauto <- sum(elast_coauto*(database$fe_pess*predictions_base_rum[,'auto']/sum(database$fe_pess*predictions_base_rum[,'auto'])))
agg_coauto <- sum(elast_coauto*(predictions_base_rum[,'auto']/sum(predictions_base_rum[,'auto'])))

#co bus
elast_cobus <- (1-predictions_base_rum[,'bus'])*model_lc$estimate[['b_co_1']]*database$co_bus
aggfv_cobus <- sum(elast_cobus*(database$fe_via*predictions_base_rum[,'bus']/sum(database$fe_via*predictions_base_rum[,'bus'])))
aggfi_cobus <- sum(elast_cobus*(database$fe_pess*predictions_base_rum[,'bus']/sum(database$fe_pess*predictions_base_rum[,'bus'])))
agg_cobus <- sum(elast_cobus*(predictions_base_rum[,'bus']/sum(predictions_base_rum[,'bus'])))

#co metro
elast_cometro <- (1-predictions_base_rum[,'metro'])*model_lc$estimate[['b_co_1']]*database$co_metro
aggfv_cometro <- sum(elast_cobus*(database$fe_via*predictions_base_rum[,'metro']/sum(database$fe_via*predictions_base_rum[,'metro'])))
aggfi_cometro <- sum(elast_cobus*(database$fe_pess*predictions_base_rum[,'metro']/sum(database$fe_pess*predictions_base_rum[,'metro'])))
agg_cometro <- sum(elast_cobus*(predictions_base_rum[,'metro']/sum(predictions_base_rum[,'metro'])))

#enumerado
agg_coauto
agg_cobus
agg_cometro
agg_ttauto
agg_ttbus
agg_ttmetro

#viagem
aggfv_coauto
aggfv_cobus
aggfv_cometro
aggfv_ttauto
aggfv_ttbus
aggfv_ttmetro

#individuo
aggfi_coauto
aggfi_cobus
aggfi_cometro
aggfi_ttauto
aggfi_ttbus
aggfi_ttmetro

### RRM
#tt auto
elast_ttauto <-   ((model_lc$estimate[['b_ttauto_2']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1) +
                   (model_lc$estimate[['b_ttauto_2']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1) +
                   ((-model_lc$estimate[['b_ttauto_2']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                   ((-model_lc$estimate[['b_ttauto_2']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                   (model_lc$estimate[['b_ttauto_2']])/(exp((model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'] +
                   (model_lc$estimate[['b_ttauto_2']])/(exp((model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$tt_auto

aggfv_ttauto <- sum(elast_ttauto*(database$fe_via*predictions_base_rrm[,'auto']/sum(database$fe_via*predictions_base_rrm[,'auto'])))
aggfv_ttauto <- sum(elast_ttauto*(database$fe_pess*predictions_base_rrm[,'auto']/sum(database$fe_pess*predictions_base_rrm[,'auto'])))
agg_ttauto <- sum(elast_ttauto*(predictions_base_rrm[,'auto']/sum(predictions_base_rrm[,'auto'])))

#tt bus
elast_ttbus <-    ((model_lc$estimate[['b_ttbus_2']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1) +
                  (model_lc$estimate[['b_ttbus_2']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1) +
                  ((model_lc$estimate[['b_ttbus_2']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                  ((-model_lc$estimate[['b_ttbus_2']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                  (model_lc$estimate[['b_ttbus_2']])/(exp((model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'auto'] +
                  (model_lc$estimate[['b_ttbus_2']])/(exp((model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$tt_bus_otp

aggfv_ttbus <- sum(elast_ttbus*(database$fe_via*predictions_base_rrm[,'bus']/sum(database$fe_via*predictions_base_rrm[,'bus'])))
aggfi_ttbus <- sum(elast_ttbus*(database$fe_pess*predictions_base_rrm[,'bus']/sum(database$fe_pess*predictions_base_rrm[,'bus'])))
agg_ttbus <- sum(elast_ttbus*(predictions_base_rrm[,'bus']/sum(predictions_base_rrm[,'bus'])))

#tt metro
elast_ttmetro <-    ((model_lc$estimate[['b_ttmetro_2']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1) +
                    (model_lc$estimate[['b_ttmetro_2']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1) +
                    ((-model_lc$estimate[['b_ttmetro_2']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'metro'] +
                    ((-model_lc$estimate[['b_ttmetro_2']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'metro'] +
                    (model_lc$estimate[['b_ttmetro_2']])/(exp((model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'auto'] +
                    (model_lc$estimate[['b_ttmetro_2']])/(exp((model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'])*database$tt_metro

aggfv_ttmetro <- sum(elast_ttmetro*(database$fe_via*predictions_base_rrm[,'metro']/sum(database$fe_via*predictions_base_rrm[,'metro'])))
aggfi_ttmetro <- sum(elast_ttmetro*(database$fe_pess*predictions_base_rrm[,'metro']/sum(database$fe_pess*predictions_base_rrm[,'metro'])))
agg_ttmetro <- sum(elast_ttmetro*(predictions_base_rrm[,'metro']/sum(predictions_base_rrm[,'metro'])))

#co auto
elast_coauto <-    ((model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1) +
                   (model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1) +
                   ((-model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                   ((-model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                   (model_lc$estimate[['b_co_2']])/(exp((model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'] +
                   (model_lc$estimate[['b_co_2']])/(exp((model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$co_auto

aggfv_coauto <- sum(elast_coauto*(database$fe_via*predictions_base_rrm[,'auto']/sum(database$fe_via*predictions_base_rrm[,'auto'])))
aggfi_coauto <- sum(elast_coauto*(database$fe_pess*predictions_base_rrm[,'auto']/sum(database$fe_pess*predictions_base_rrm[,'auto'])))
agg_coauto <- sum(elast_coauto*(predictions_base_rrm[,'auto']/sum(predictions_base_rrm[,'auto'])))

#co bus
elast_cobus <-    ((model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1) +
                  (model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1) +
                  ((-model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                  ((-model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                  (model_lc$estimate[['b_co_2']])/(exp((model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$ estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'auto'] +
                  (model_lc$estimate[['b_co_2']])/(exp((model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$co_bus

aggfv_cobus <- sum(elast_cobus*(database$fe_via*predictions_base_rrm[,'bus']/sum(database$fe_via*predictions_base_rrm[,'bus'])))
aggfi_cobus <- sum(elast_cobus*(database$fe_pess*predictions_base_rrm[,'bus']/sum(database$fe_pess*predictions_base_rrm[,'bus'])))
agg_cobus <- sum(elast_cobus*(predictions_base_rrm[,'bus']/sum(predictions_base_rrm[,'bus'])))

#co metro
elast_cometro <-    ((model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1) +
                    (model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1) +
                    ((-model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'metro'] +
                    ((-model_lc$estimate[['b_co_2']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'metro'] +
                    (model_lc$estimate[['b_co_2']])/(exp((model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'auto'] +
                    (model_lc$estimate[['b_co_2']])/(exp((model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'])*database$co_metro

aggfv_cometro <- sum(elast_cobus*(database$fe_via*predictions_base_rrm[,'metro']/sum(database$fe_via*predictions_base_rrm[,'metro'])))
aggfi_cometro <- sum(elast_cobus*(database$fe_pess*predictions_base_rrm[,'metro']/sum(database$fe_pess*predictions_base_rrm[,'metro'])))
agg_cometro <- sum(elast_cobus*(predictions_base_rrm[,'metro']/sum(predictions_base_rrm[,'metro'])))

#enumerado
agg_coauto
agg_cobus
agg_cometro
agg_ttauto
agg_ttbus
agg_ttmetro

#viagem
aggfv_coauto
aggfv_cobus
aggfv_cometro
aggfv_ttauto
aggfv_ttbus
aggfv_ttmetro

#individuo
aggfi_coauto
aggfi_cobus
aggfi_cometro
aggfi_ttauto
aggfi_ttbus
aggfi_ttmetro

### RRM PAPER
#tt auto
elast_ttauto <- ((model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1) +
                   (model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1) +
                   ((model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                   ((model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                   (model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'] +
                   (model_lc$estimate[['b_ttauto_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttauto_2']]*database$tt_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$tt_auto

aggfv_ttauto <- sum(elast_ttauto*(database$fe_via*predictions_base_rrm[,'auto']/sum(database$fe_via*predictions_base_rrm[,'auto'])))
aggfv_ttauto <- sum(elast_ttauto*(database$fe_pess*predictions_base_rrm[,'auto']/sum(database$fe_pess*predictions_base_rrm[,'auto'])))
agg_ttauto <- sum(elast_ttauto*(predictions_base_rrm[,'auto']/sum(predictions_base_rrm[,'auto'])))

#tt bus
elast_ttbus <- ((model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1) +
                  (model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1) +
                  ((model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                  ((model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                  (model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'] +
                  (model_lc$estimate[['b_ttbus_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttmetro_2']]*database$tt_metro - model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$tt_bus_otp

aggfv_ttbus <- sum(elast_ttbus*(database$fe_via*predictions_base_rrm[,'bus']/sum(database$fe_via*predictions_base_rrm[,'bus'])))
aggfi_ttbus <- sum(elast_ttbus*(database$fe_pess*predictions_base_rrm[,'bus']/sum(database$fe_pess*predictions_base_rrm[,'bus'])))
agg_ttbus <- sum(elast_ttbus*(predictions_base_rrm[,'bus']/sum(predictions_base_rrm[,'bus'])))

#tt metro
elast_ttmetro <- ((model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1) +
                    (model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1) +
                    ((model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                    ((model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                    (model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttauto_2']]*database$tt_auto - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'] +
                    (model_lc$estimate[['b_ttmetro_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_ttbus_2']]*database$tt_bus_otp - model_lc$estimate[['b_ttmetro_2']]*database$tt_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'])*database$tt_metro

aggfv_ttmetro <- sum(elast_ttmetro*(database$fe_via*predictions_base_rrm[,'metro']/sum(database$fe_via*predictions_base_rrm[,'metro'])))
aggfi_ttmetro <- sum(elast_ttmetro*(database$fe_pess*predictions_base_rrm[,'metro']/sum(database$fe_pess*predictions_base_rrm[,'metro'])))
agg_ttmetro <- sum(elast_ttmetro*(predictions_base_rrm[,'metro']/sum(predictions_base_rrm[,'metro'])))

#co auto
elast_coauto <- ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1) +
                   (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1) +
                   ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                   ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'bus'] +
                   (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'] +
                   (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_auto)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$co_auto

aggfv_coauto <- sum(elast_coauto*(database$fe_via*predictions_base_rrm[,'auto']/sum(database$fe_via*predictions_base_rrm[,'auto'])))
aggfi_coauto <- sum(elast_coauto*(database$fe_pess*predictions_base_rrm[,'auto']/sum(database$fe_pess*predictions_base_rrm[,'auto'])))
agg_coauto <- sum(elast_coauto*(predictions_base_rrm[,'auto']/sum(predictions_base_rrm[,'auto'])))

#co bus
elast_cobus <- ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1) +
                  (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1) +
                  ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                  ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                  (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$ estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'] +
                  (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_metro - model_lc$estimate[['b_co_2']]*database$co_bus)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'metro'])*database$co_bus

aggfv_cobus <- sum(elast_cobus*(database$fe_via*predictions_base_rrm[,'bus']/sum(database$fe_via*predictions_base_rrm[,'bus'])))
aggfi_cobus <- sum(elast_cobus*(database$fe_pess*predictions_base_rrm[,'bus']/sum(database$fe_pess*predictions_base_rrm[,'bus'])))
agg_cobus <- sum(elast_cobus*(predictions_base_rrm[,'bus']/sum(predictions_base_rrm[,'bus'])))

#co metro
elast_cometro <- ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1) +
                    (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp(-(model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1) +
                    ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                    ((model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1))*predictions_base_rrm[,'auto'] +
                    (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_auto - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'] +
                    (model_lc$estimate[['b_co_2']]/model_lc$estimate[['b_mu']])/(exp((model_lc$estimate[['b_co_2']]*database$co_bus - model_lc$estimate[['b_co_2']]*database$co_metro)/model_lc$estimate[['b_mu']])+1)*predictions_base_rrm[,'bus'])*database$co_metro

aggfv_cometro <- sum(elast_cobus*(database$fe_via*predictions_base_rrm[,'metro']/sum(database$fe_via*predictions_base_rrm[,'metro'])))
aggfi_cometro <- sum(elast_cobus*(database$fe_pess*predictions_base_rrm[,'metro']/sum(database$fe_pess*predictions_base_rrm[,'metro'])))
agg_cometro <- sum(elast_cobus*(predictions_base_rrm[,'metro']/sum(predictions_base_rrm[,'metro'])))

#enumerado
agg_coauto
agg_cobus
agg_cometro
agg_ttauto
agg_ttbus
agg_ttmetro

#viagem
aggfv_coauto
aggfv_cobus
aggfv_cometro
aggfv_ttauto
aggfv_ttbus
aggfv_ttmetro

#individuo
aggfi_coauto
aggfi_cobus
aggfi_cometro
aggfi_ttauto
aggfi_ttbus
aggfi_ttmetro
 

