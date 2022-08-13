# packages
library(tidyverse)
library(nnet)
library(caret)
library(smd)
library(kableExtra)
# fonction
source("Source_papier_prett.R",
       encoding = "UTF-8"
)

# var
Id <- paste0(quote(id_patient)) # change from recently if bug here "id_patient"
name_expo <- paste0(quote(INS_obs_categ))
name_outcome <- paste0(quote(SYM_echellePichot))
visite_seuil <- 2
# data
df_modele_pds <- read_rds("data/genere/data_pre_impute_Pichot_nb_vis_2.rds") 



# perform model variables selection

select_var <- read_rds("data/genere/model_outcome.rds")


remove_var_model <- c(
  "SYM_somnolenceDiurne",
  "Vis_init_PSG_saturationMoyenne","Vis_init_PSG_delai_explo"
  ,"PPC_pressionMoyenne","INS_PPC_effet_ind","PPC_id_typeMasque"
)


Var_model_outcome <- select_var$var_select$model_outcome[select_var$var_select$model_outcome %notin% remove_var_model]





# function that compute outcome estimation with linear regression and IPWRA
# and compute weight  IPTW  for IPW
function_res_papier <- function(data, index, params_fun) {
  df_fun <- data[index, ]
  
  Id_patient_analyse <- df_fun %>% select(all_of(c(params_fun$Id_data, params_fun$expo)))
  
  res_function_poids <- do.call(calcule_pds_stage, c(
    list(donne = df_fun),
    params_fun
  ))
  
  
  
  poids_strap <- res_function_poids$res_intermediaire$poids$poids_tronc$poids_trunc %>%
    unlist(use.names = FALSE)
  
  
  donnee_pre_mod_fin <- df_fun %>%
    mutate_at(params_fun$expo, list(~ relevel(., ref = "4"))) %>%
    select(-all_of(params_fun$Id_data))
  
  model_pred_outcomes <- glm(as.formula(paste0(params$out_come, " ~ ",
                                               paste0(params_fun$expo," + "),
                                               paste0(Var_model_outcome,collapse = " + ")
  )),
  data = donnee_pre_mod_fin,
  family = gaussian(),
  weights = poids_strap
  )
  
  
  model_pred_outcomes_non_pond <- glm(as.formula(paste0(params$out_come, " ~ ",
                                                        paste0(params_fun$expo," + "),
                                                        paste0(Var_model_outcome,collapse = " + ")
  )),
  data = donnee_pre_mod_fin,
  family = gaussian()
  )
  IPTW_effect <- glm(as.formula(paste0(params$out_come, " ~ ",params_fun$expo
  )),
  data = donnee_pre_mod_fin,
  family = gaussian(),
  weights = poids_strap
  )
  
  df_pred_ss_pond <- predict(model_pred_outcomes_non_pond, newdata = donnee_pre_mod_fin)
  
  unique_expo <- df_fun %>%
    select(all_of(params_fun$expo)) %>%
    unlist(use.names = FALSE)
  
  expo_fact_fun <- params_fun$expo
  
  df_pred <- donnee_pre_mod_fin %>%
    select(-all_of(params_fun$expo)) %>% # data.table::
    merge(data.frame(temp_name = unique(unique_expo))) %>%
    rename(!!expo_fact_fun := temp_name)
  
  treatment_effect_df <- predict(model_pred_outcomes, newdata = df_pred) %>% matrix(nrow = nrow(df_fun), ncol = length(unique(unique_expo)))
  
  colnames(treatment_effect_df) <- paste0("treatment_grp_", unique(unique_expo))
  
  return(list(
    var_selec = Var_model_outcome,
    poids = poids_strap,
    reg_pond = coef(model_pred_outcomes),
    Treatment_effect = treatment_effect_df,
    df_id_grp_tt = Id_patient_analyse,
    result_func = res_function_poids$df,
    reg_non_pond = coef(model_pred_outcomes_non_pond),
    estim_reg_lin = df_pred_ss_pond,
    ind = index,
    mod_poids_coef = coef(res_function_poids$res_intermediaire$regression_modele_pds$regression_temps_ind),
    IPTW_res = coef(IPTW_effect)
  ))
}


params <- list(
  expo = name_expo,
  covar = Var_model_outcome, # var_mod_pds ,
  Id_data = Id,
  weighting_function = multinomial_IPTW,
  out_come = name_outcome,
  percentile_tronc = 0.01
)
donnee_boot <- select_var$data_frame$df_analyse %>%
  mutate_at(
    typages_function(select_var$data_frame$df_analyse, 10)$colonne_type$facteur,
    factor
  )


if ("OK" == "pas OK") {
  boot_strap_full_process <- boot_strap_fun(donnee_boot, function_res_papier, 1000, params)
  saveRDS(boot_strap_full_process, file = paste0("data/genere/bott_strap_full_process_pichot.rds"))
}



boot_strap_full_process <- readRDS("data/genere/bott_strap_full_process_pichot.rds")


###############################################################################
# boot strap analysis
# truncated weight
boot_strap_pds_trunc <- lapply(boot_strap_full_process[6, ], function(x) {
  Trucature_pds_function(x$Weight, c(0, 1, 5, 10, 25, 50) / 100) %>%
    summarise_all(list(
      ~ mean(.),
      ~ min(.),
      ~ max(.)
    )) %>%
    t.df() %>%
    separate(key,
             into = c(
               "truncations",
               "fun"
             ),
             sep = "_"
    ) %>%
    pivot_wider(
      names_from = fun,
      values_from = col_1
    ) %>%
    rename_all(function(x) c("Truncations", "Mean", "Minimum", "Maximum"))
}) %>%
  bind_rows()
summarry_boot_strap_pds_trunc <- boot_strap_pds_trunc %>%
  mutate(Truncations = factor(Truncations,
                              levels = c(
                                "(0; 1)",
                                "(0.01; 0.99)",
                                "(0.05; 0.95)",
                                "(0.1; 0.9)",
                                "(0.25; 0.75)",
                                "(0.5; 0.5)"
                              )
  )) %>%
  group_by(Truncations) %>%
  summarise_at(vars(-group_cols()), ~ paste0(
    sprintf("%.01f", round(mean(.), 1)),
    "(", sprintf("%.01f", round(quantile(., 0.05), 1)), "; ",
    sprintf("%.01f", round(quantile(., 0.95), 1)), ")"
  ))

# Weight summary among bootstrap
descriptif_pds_boot_sumary_table <- lapply(boot_strap_full_process[6, ], function(x) {
  x %>%
    group_by(exposition) %>%
    summarise(
      "mean_W" = mean(Weight),
      "min_W" = min(Weight),
      "max_W" = max(Weight),
      "mean_SW" = mean(SWeight),
      "min_SW" = min(SWeight),
      "max_SW" = max(SWeight),
      .groups = "drop"
    )
}) %>% bind_rows()


summarydescriptif_pds_boot <- descriptif_pds_boot_sumary_table %>%
  group_by(exposition) %>%
  # mutate_all(list(~as.character())) %>%
  summarise_at(vars(-group_cols()), ~ paste0(sprintf("%.01f", round(mean(.), 1)), "(", sprintf("%.01f", round(quantile(., 0.025), 1)), "; ", sprintf("%.01f", round(quantile(., 0.975), 1)), ")"))


# Bootstrap comparaison of treatment effect estimation

boot_strap_coef_reg_non_pond <- apply(boot_strap_full_process, 2, function(x) {
  weight_fun <- x[2]$poids
  id_boot_strap <- x[5]$df_id_grp_tt %>% select(-INS_obs_categ)
  df_diff_simple <- id_boot_strap %>% left_join(donnee_boot, by = c("id_patient")) # %>% select(-id_patient)
  rbind(
    x[7]$reg_non_pond %>% as.data.frame() %>% rename("value" = 1) %>%
      rownames_to_column() %>% filter(str_detect(rowname, "INS_obs_categ")) %>% 
      t.df(pivot = "rowname") %>%
      mutate(key = "reg_pond"),
    x[11]$IPTW_res %>% as.data.frame() %>% rename("value" = 1) %>%
      rownames_to_column() %>% filter(str_detect(rowname, "INS_obs_categ")) %>% 
      t.df(pivot = "rowname") %>%
      mutate(key = "IPTW"),
    x[3]$reg_pond %>% as.data.frame() %>%
      rename("value" = 1) %>% rownames_to_column() %>%
      filter(str_detect(rowname, "INS_obs_categ")) %>%
      t.df(pivot = "rowname") %>% mutate(key = "reg_non_pond"),
    df_diff_simple %>% mutate(weight = weight_fun) %>%
      group_by(INS_obs_categ) %>%
      summarise(
        normal_mean = mean(SYM_echellePichot),
        # IPTW = weighted.mean(x = SYM_echelleEpworth, w = weight)
      ) %>%
      t.df("INS_obs_categ") %>%
      mutate_at(vars(-key), ~ (. - `4`)) %>%
      select(-`4`) %>%
      rename_at(vars(-key), ~ paste0("INS_obs_categ", .))
  )
}) %>% bind_rows()


for (i in 1:3) {
  print(c(
    mean(boot_strap_coef_reg_non_pond %>% filter(key == "IPTW") %>%
           select(contains(as.character(i))) %>% unlist(use.names = FALSE)) -
      mean(boot_strap_coef_reg_non_pond %>% filter(key == "reg_pond") %>%
             select(contains(as.character(i))) %>% unlist(use.names = FALSE))
  ))
}

t_test_diff_ess_bygroup_iptw_piwra <- sapply(1:3, function(i) {
  t.test(
    boot_strap_coef_reg_non_pond %>% filter(key == "IPTW") %>%
      select(contains(as.character(i))) %>% unlist(use.names = FALSE),
    boot_strap_coef_reg_non_pond %>% filter(key == "reg_pond") %>%
      select(contains(as.character(i))) %>% unlist(use.names = FALSE)
  )
})






df_4meth <- boot_strap_coef_reg_non_pond %>%
  pivot_longer(-key) %>%
  group_by(key, name) %>%
  summarise(moy = mean(value), ci_sup = quantile(value, 0.975),
            ci_inf = quantile(value, 0.025),
            .groups = "drop") %>%
  mutate(
    key = str_replace(key, "normal_mean", "Unadjusted mean comparaison"),
    key = str_replace(key, "reg_non_pond", "Multivariable linear regression"),
    key = str_replace(key, "reg_pond", "IPWRA")
  ) %>%
  mutate(name = str_remove(name, "INS_obs_categ")) %>%
  mutate(key = factor(key, levels = c("Unadjusted mean comparaison", "IPTW", "Multivariable linear regression", "IPWRA")))


df_4meth %>%
  #filter(key == "IPTW" | key == "IPWRA") %>%
  mutate(ci = paste0("(", round(ci_inf, 2), "; ", round(ci_sup, 2), ")"), moy = round(moy, 2)) %>%
  select(-ci_sup, -ci_inf) %>%
  pivot_wider(names_from = key, values_from = c(moy, ci)) 
##############################################################################



list_id_boot <- lapply(boot_strap_full_process[5, ], function(x) {
  x$id_patient[x$id_patient %notin% donnee_boot$id_patient]
})



##############################################################################
# Bootstrap Variables balance before and after weighting
nb_grp <- donnee_boot %>%
  select(all_of(name_expo)) %>%
  unique() %>%
  unlist(use.names = FALSE) %>%
  sort()


df_smd <- select_var$data_frame$df_analyse %>% select(all_of(c(Id, name_expo, Var_model_outcome))) # one_hot_fb(select_var$data_frame$df_analyse ,list_factor =  c(  "PPC_id_typeMasque"  ))
descriptif_pds_boot <- apply(boot_strap_full_process, 2, function(i) {
  weight_fun <- i[2]$poids
  id_boot_strap <- i[5]$df_id_grp_tt %>% select(-INS_obs_categ)
  
  
  
  df_func <- id_boot_strap %>% left_join(df_smd, by = c("id_patient")) # %>% select(-id_patient)
  
  mat_var <- df_func %>%
    select(-id_patient, -all_of(c(name_expo))) %>%
    as.matrix()
  
  
  
  df_max_func <- data.frame(
    sm_max_non_pondere = apply(smd(mat_var,
                                   df_func$INS_obs_categ),
                               2, function(x) abs(mean(x))
    ),
    sm_max_pondere = apply(smd(mat_var, df_func$INS_obs_categ, weight_fun), 2, function(x) abs(mean(x)))
  ) #
})

df_max <- apply(simplify2array(lapply(descriptif_pds_boot, as.matrix)), 1:2, mean) %>%
  as.data.frame() %>%
  rownames_to_column()


#################################################################################
IPWRA_coef_boot <- lapply(boot_strap_full_process[3, ], function(x) {
  x
}) %>% bind_rows()

IPWRA_coef_boot_CI <- IPWRA_coef_boot %>%
  select(-`(Intercept)`) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    moy = mean(value),
    ci_sup = quantile(value, 0.975),
    ci_inf = quantile(value, 0.025),
    .groups = "drop"
  )

#################################################################################
weight_coef_boot <- lapply(boot_strap_full_process[10, ], function(x) {
  as.data.frame(x)
}) %>% bind_rows()

weight_coef_boot_CI <- weight_coef_boot %>% 
  rownames_to_column() %>%
  select(-`(Intercept)`) %>%
  mutate(rowname = str_extract(rowname, "^\\d{1}")) %>%
  pivot_longer(-rowname) %>%
  group_by(name, rowname) %>%
  summarise(
    moy = mean(value),
    ci_sup = quantile(value, 0.975),
    ci_inf = quantile(value, 0.025),
    .groups = "drop"
  )

#################################################################################
reg_coef_boot <- lapply(boot_strap_full_process[7, ], function(x) {
  as.data.frame(x) %>% rownames_to_column() %>%  t.df(pivot="rowname") 
})  %>% bind_rows() %>% select(-key)

reg_coef_boot_CI <- reg_coef_boot %>%
  select(any_of(colnames(weight_coef_boot)),contains("INS_obs_categ")) %>% 
  select(-`(Intercept)`) %>%
  pivot_longer(everything()) %>%
  group_by(name) %>%
  summarise(
    moy = mean(value),
    ci_sup = quantile(value, 0.975),
    ci_inf = quantile(value, 0.025),
    .groups = "drop"
  )





##################################################################################
# Print paper result
vec_cut <- read_rds("data/genere/obs_cut_pichot2.rds")

# patient per group and our range
for (i in seq_along(vec_cut[-1])) {
  cat("\\item", table(df_modele_pds$INS_obs_categ)[i],
      " ( ",
      round(
        table(df_modele_pds$INS_obs_categ)[i] * 100 / sum(table(df_modele_pds$INS_obs_categ)),
        1), "\\% )", " patients with an adherence between ",
      roud_hour(vec_cut)[i], " and ", roud_hour(vec_cut)[i + 1], "\n")
}

groupe_obs_table_latex <- sapply(seq_along(vec_cut[-1]), function(i) (paste0(vec_cut[i], "-", vec_cut[i + 1], " h")))




alpha_tableau_resume <- 0.05
var_instrumental_name <- c(
  "INS_obs_categ" = "adherence groups",
  "Vis_init_INS_IAH" = "apnea hypopnea index",
  "INS_PPC_effet_ind" = "number of ADR types under CPAP",
  "INS_tmp_entre_rdv" = "duration since diagnosis (year)"
)




table_rename <- one_hot_fb(select_var$data_frame$df_analyse, list_factor = c("PPC_id_typeMasque")) %>%
  select(-Id) %>%
  rename_variables(var_instrum_name = var_instrumental_name)



# nb var model fin


print_var_at_diag <- str_replace(
  str_replace(
    str_replace(
      str_remove_all(
        str_replace(
          str_replace(
            str_replace(
              firstlow(table_rename$complete_table$Label[table_rename$complete_table$var_name %in% Var_model_outcome[str_detect(
                Var_model_outcome,
                "^Vis_init_|^Mesure_unique_"
              )]]),
              "epworth sleepiness scale", "daytime sleepiness measured by ESS scale"
            ),
            "pichot's fatigue scale", "fatigue measured by Pichot's scale"
          ),
          "depression scale", "depression  measured by Pichot's depression scale"
        ),
        "diagnosis" #|\\(.+\\)"
      ), "\\(kg/m2\\)", "(kg/m\textsuperscript{2})"
    ),
    "neck circumference", "neck circumference (cm)"
  ), # add into variables dictionnary instead
  "SaO2", "SaO2 (\\\\%)"
) # add into variables dictionnary instead


cat(str_replace(paste0("are the following variables at diagnosis : ", paste0(print_var_at_diag, collapse = ", ")), ",(?=[^,]*$)", " and"))

cat(str_replace(paste0("The variables under CPAP treatment : ", paste0(str_replace(
  str_replace(
    str_replace(
      firstlow(table_rename$complete_table$Label[table_rename$complete_table$var_name %in% Var_model_outcome[!str_detect(Var_model_outcome, "^Vis_init_|^Mesure_unique_")]]),
      "pichot's fatigue scale", "fatigue measured by Pichot's scale"
    ),
    "ADR", "adverse drug reaction"
  ),
  "depression scale", "depression  measured by Pichot's depression scale"
),
collapse = ", "
)), ",(?=[^,]*$)", " and"))



table_rename$complete_table$var_name[table_rename$complete_table$var_name %in% Var_model_outcome[!str_detect(Var_model_outcome, "^Vis_init_|^Mesure_unique_")]]


select_var$data_frame$df_analyse %>%
  summarise(
    epworth = mean(SYM_echellePichot),
    sd = sd(SYM_echellePichot),
    base_epworth = mean(Vis_init_SYM_echellePichot),
    Sd_base = sd(Vis_init_SYM_echellePichot),
    diff_epworth = median(Vis_init_SYM_echellePichot - SYM_echellePichot),
    SD_diff = sd(Vis_init_SYM_echellePichot - SYM_echellePichot)
  ) %>%
  arrondie_df(2)


select_var$data_frame$df_analyse %>%
  group_by(INS_obs_categ) %>%
  summarise(
    epworth = mean(SYM_echellePichot),
    SD = sd(SYM_echellePichot),
    base_epworth = mean(Vis_init_SYM_echellePichot),
    Sd_base = sd(Vis_init_SYM_echellePichot),
    diff_epworth = mean(Vis_init_SYM_echellePichot - SYM_echellePichot),
    SD_diff = sd(Vis_init_SYM_echellePichot - SYM_echellePichot),
  ) %>%
  arrondie_df(2)


quantile(select_var$data_frame$df_analyse$Vis_init_SYM_echellePichot, c(0.25, 0.75))


# number of patient with severe osa
sum(select_var$data_frame$df_analyse$Vis_init_INS_IAH >= 30)
round(mean(select_var$data_frame$df_analyse$Vis_init_INS_IAH >= 30) * 100, 1)


# Final result
df_4meth %>%
  mutate(ci = paste0("(", round(ci_inf, 2), "; ", round(ci_sup, 2), ")"), moy = round(moy, 2)) %>%
  select(-ci_sup, -ci_inf) %>%
  pivot_wider(names_from = key, values_from = c(moy, ci))



# Anova of initial epworth scale
Init_ESS_mean_comp <- aov(Vis_init_SYM_echellePichot ~ INS_obs_categ, donnee_boot)
summary(Init_ESS_mean_comp)


# Anova simple mean follow-up Epworth
Folow_ESS_mean_comp <- aov(SYM_echellePichot ~ INS_obs_categ, donnee_boot)
summary(Folow_ESS_mean_comp)

# mean initial obs for all patients
obs_nume <- read_rds("data/genere/obs_pichot_quantit2.rds")
obs_nume %>%
  inner_join(donnee_boot, by = Id) %>%
  summarise(roud_hour(mean(PPC_observanceMoy_finale)))


# pval diff obs pond


ESS_test_by_group_AIPW <- boot_strap_coef_reg_non_pond %>%
  filter(key == "reg_pond")
ESS_test_by_group_IPWT <- boot_strap_coef_reg_non_pond %>%
  filter(key == "IPTW")
table_col <- combn(ESS_test_by_group_AIPW %>% select(-key) %>% colnames(), 2)


mini_test_fun <- function(df, table_col) {
  unslit_fb <- function(Y) unlist(Y, use.names = FALSE)
  
  apply(table_col, 2, function(x) {
    t.test(unslit_fb(df[x[1]]), unslit_fb(df[x[2]]))
  })
}





mini_test_fun(ESS_test_by_group_AIPW, table_col)


mini_test_fun(ESS_test_by_group_IPWT, table_col)

# 
# 
# ##################################################################################
# # Figure and table generation
# 
# # Netoyage du fichier tex
# path_latex <- "figure_papier/table_iptw.tex"
# 
# sub_file_latex_preambule <- {
#   "\\documentclass[../main.tex]{subfiles}
# \\begin{document}"
# }
# 
# 
# 
# write(sub_file_latex_preambule, path_latex, append = FALSE)
# # génération du code pour le papier inutile dans le rapport
# summarydescriptif_pds_boot %>%
#   select(exposition, ends_with("_W")) %>%
#   rename(
#     `Adherence group` = exposition,
#     Mean = mean_W,
#     Minimum = min_W,
#     Maximum = max_W
#   ) %>%
#   mutate(`Adherence group` = groupe_obs_table_latex) %>%
#   kable(
#     format = "latex",
#     align = "c",
#     booktabs = TRUE,
#     escape = FALSE,
#     label = "Unstabilized_table",
#     linesep = "",
#     caption = "Distribution of weights"
#   ) %>%
#   kable_styling(latex_options = c("striped", "HOLD_position", "scale_down", "repeat_header")) %>%
#   footnote(general = c("Data are presented as mean (95% confidence interval) of bootstrap iterations"), general_title = "") %>%
#   ecriture_fichier(path_latex)
# write("\\clearpage", path_latex, append = TRUE)
# 
# 
# 
# # in paper we reduct number of variables in tables
# Variables_core_paper <- c(
#   "INS_obs_categ", "SYM_echelleEpworth", "Vis_init_SYM_echelleEpworth", "Mesure_unique_CHA_id_sexe", "Vis_init_CHA_age", "Vis_init_CHA_BMI", "Vis_init_INS_IAH", "PPC_IRAH", "Vis_init_FDR_fumeur",
#   "Vis_init_SYM_echelleDepression", "SYM_echelleDepression", "Vis_init_SYM_echellePichot", "SYM_echellePichot",
#   "Vis_init_SYM_cephaleesMatinales", "SYM_cephaleesMatinales", "Vis_init_SYM_fatigueMatinale", "SYM_fatigueMatinale",
#   "Mesure_unique_ATC_diabete", "INS_PPC_effet_ind", "INS_tmp_entre_rdv"
# )
# foot_note <- list(
#   CPAP = list(var = c("INS_PPC_effet_ind", "PPC_IRAH"), foot_notes = "CPAP : Continuous Positive Airway Pressure"),
#   ADR = list(var = c("INS_PPC_effet_ind"), foot_notes = "ADR : adverse drug reaction")
# )
# 
# 
# table__resume_rapport_non_rename <- select_var$data_frame$df_analyse %>%
#   select(any_of(Variables_core_paper))
# 
# 
# foot_note_var_presence <- lapply(foot_note, function(x) {
#   x$foot_notes[any(x$var %in% colnames(table__resume_rapport_non_rename))]
# }) %>% unlist()
# # Tables avec les variables renommer pour jointure
# 
# 
# table_resume_rapport_rename <- table__resume_rapport_non_rename %>% rename_variables(var_instrum_name = var_instrumental_name)
# table__resume_rapport <- table_resume_rapport_rename$table_rename
# 
# 
# diag_var_table_latex <- table_resume_rapport_rename$complete_table %>%
#   mutate(diagnostic = !is.na(preffix_all)) %>%
#   select(Label, diagnostic)
# 
# 
# table_resume_rapport_latex_ss_escape <- table_resume_latex(
#   df_one_hot_encode_fonc = table__resume_rapport,
#   name_expo_fonc = "Adherence groups",
#   nom_grp = "Adherence grp",
#   p_val = TRUE,
#   alpha = alpha_tableau_resume
# )
# 
# table_resume_rapport_latex_ss_escape_order <- table_resume_rapport_latex_ss_escape %>%
#   rename_all(~ c(unlist(table_resume_rapport_latex_ss_escape["Number of patient", ], use.names = FALSE))) %>%
#   # Marquage modèle de poids
#   rownames_to_column() %>%
#   filter(rowname != "Number of patient") %>%
#   inner_join(diag_var_table_latex, by = c("rowname" = "Label")) %>%
#   mutate(rowname = R.utils::capitalize(str_remove(rowname, "^Diagnosis "))) %>%
#   arrange(desc(diagnostic))
# 
# 
# write("\\clearpage", path_latex, append = TRUE)
# header_landscape <- "\\newgeometry{a4paper,left=1in,right=1in,top=1in,bottom=1in,nohead}" # Il faut changer les marges des pages avant de landscape pour éviter un saut de pages indésiré
# write(header_landscape, path_latex, append = TRUE)
# 
# 
# # remove_rownames %>%
# # column_to_rownames(var = "rowname") %>%
# table_resume_rapport_latex_ss_escape_order %>%
#   select(-diagnostic) %>%
#   mutate(
#     rowname = str_replace(rowname, "Epworth sleepiness scale", "ESS score"),
#     rowname = str_replace(rowname, "Gender", "Gender (male)"),
#     rowname = str_replace(rowname, "Apnea hypopnea index", " Apnea hypopnea index (event/h)")
#   ) %>%
#   kable(
#     format = "latex",
#     align = "c",
#     booktabs = TRUE,
#     escape = FALSE,
#     label = "Table_resume_pop",
#     linesep = "",
#     col.names = c("", colnames(table_resume_rapport_latex_ss_escape_order)[c(-1, -ncol(table_resume_rapport_latex_ss_escape_order))]),
#     caption = "Patient characteristics according to the adherence group"
#   ) %>%
#   add_header_above(c("", "All groups", paste0(groupe_obs_table_latex, " (", seq_along(groupe_obs_table_latex), ")")),
#                    line = FALSE
#   ) %>%
#   pack_rows(
#     "Variables at diagnosis ", 1,
#     sum(table_resume_rapport_latex_ss_escape_order$diagnostic)
#   ) %>%
#   pack_rows(
#     "Variables at follow-up",
#     sum(table_resume_rapport_latex_ss_escape_order$diagnostic) + 1,
#     nrow(table_resume_rapport_latex_ss_escape_order)
#   ) %>%
#   kable_styling(latex_options = c(
#     "striped", "HOLD_position",
#     "repeat_header"
#   ), ) %>%
#   footnote(
#     general = c(paste0(c("ESS : Epworth Sleepiness Scale", foot_note_var_presence), collapse = "; ")), general_title = "",
#     symbol = c(
#       "Quantitative variables are presented as mean (standard deviation). Qualitative variables are expressed in number of individuals (% of individuals) ",
#       "t-test was performed for the quantitative variables and a Pearson's Chi squared test for the categorical variables after application of a Bonferroni correction for multiple testing.",
#       "1,2,3,4 numbers in subscript refers to columns statistically different at the 5% threshold. e.g. 1 means that there is a statistically significant difference between the  adherence group \\textcolor{red}{of that column} and the 0-4h adherence group (1) for the variable in question."
#       # "* after index corresponding to variables with too few people to perform test"
#     ),
#     threeparttable = TRUE
#   ) %>%
#   landscape() %>%
#   ecriture_fichier(path_latex)
# end_landscape <- "\\restoregeometry % Restore the global document page margins"
# write(end_landscape, path_latex, append = TRUE)
# write("\\clearpage", path_latex, append = TRUE)
# 
# 
# 
# 
# 
# df_box_plot_fig2 <- select_var$data_frame$df_analyse %>%
#   rename_variables(var_instrum_name = var_instrumental_name) %>%
#   .$table_rename %>%
#   select(`Adherence groups`, `Diagnosis epworth sleepiness scale`, `Epworth sleepiness scale`) %>%
#   rename("Epworth sleepiness scale at diagnosis" = "Diagnosis epworth sleepiness scale", "Epworth sleepiness scale at follow-up visit" = "Epworth sleepiness scale") %>%
#   pivot_longer(-`Adherence groups`) %>%
#   mutate(
#     `Adherence groups` = as.factor(as.character(`Adherence groups`)),
#     name = factor(name, levels = c(
#       "Epworth sleepiness scale at diagnosis",
#       "Epworth sleepiness scale at follow-up visit"
#     ))
#   )
# 
# 
# ggplot(
#   data = df_box_plot_fig2,
#   aes(x = name, y = value, fill = `Adherence groups`)
# ) +
#   geom_boxplot() +
#   ggsci::scale_fill_lancet(labels = groupe_obs_table_latex) +
#   theme(legend.key = element_blank()) +
#   labs(
#     y = "Epworth sleepiness score",
#     x = paste0(""),
#     fill = "CPAP adherence"
#   )
# ggsave("figure_papier/figure_3.jpeg",
#        width = 8,
#        height = 8,
#        units = c("in")
# )
# 
# 
# # SD bootstrap
# 
# df_smd_plot <- df_max %>%
#   filter(rowname %in% Var_model_outcome)
# # filter(rowname %in% Variables_core_paper)
# df_plot_max <- table_rename$complete_table %>%
#   select(var_name, Label) %>%
#   right_join(df_smd_plot, by = c("var_name" = "rowname")) %>%
#   select(-var_name) %>%
#   mutate(Label = str_replace(
#     Label, "Diagnosis apnea hypopnea index",
#     "Apnea Hypopnea Index at diagnosis"
#   )) %>%
#   mutate(Label = factor(Label,
#                         levels = .$Label[order(abs(df_smd_plot$sm_max_non_pondere))]
#   )) %>%
#   pivot_longer(-Label) %>%
#   mutate(
#     name = str_replace_all(name, "sm_max_non_pondere", "SMD without Weighting"),
#     `Weighting` = str_replace_all(name, "sm_max_pondere", "SMD with Weighting")
#   ) %>%
#   select(-name)
# 
# ggplot(
#   data = df_plot_max,
#   mapping = aes(x = Label, y = value, group = Weighting, color = Weighting)
# ) +
#   geom_point() +
#   geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
#   scale_y_continuous(breaks = sort(c(seq(0, 1, length.out = 5), 0.1))) +
#   coord_flip() +
#   theme_bw() +
#   theme(
#     legend.background = element_rect(fill = "transparent", colour = "transparent"),
#     legend.key = element_blank(),
#     legend.title = element_blank(),
#     legend.position = c(.75, 0.15)
#   ) +
#   ggplot2::labs(
#     title = "", y = "SMD",
#     x = ""
#   )
# 
# ggsave("figure_papier/figure_4.jpeg",
#        width = 8,
#        height = 8,
#        units = c("in")
# )
# 
# 
# 
# pd <- position_dodge(0.1) # move them .05 to the left and right
# 
# df_figure_5 <- df_4meth %>%
#   #filter(key != "Unweighted linear regression") %>%
#   mutate(name = groupe_obs_table_latex[as.numeric(name)]) %>%
#   mutate(key = recode_factor(key, "Unadjusted mean comparaison" = "Unadjusted"))
# 
# 
# 
# 
# ggplot(
#   df_figure_5,
#   aes(x = name, y = moy)
# ) +
#   geom_errorbar(aes(ymin = ci_inf, ymax = ci_sup), colour = "black", width = .1, position = pd) +
#   geom_point(position = pd, size = 3) +
#   facet_wrap(~key) +
#   labs(x = "Adherence groups", y = "Change in Epworth score") +
#   theme(
#     strip.text = element_text(
#       size = 12,
#       face = "bold"
#     ) # ,
#     # text = element_text(size = 20)
#   )
# 
# 
# 
# 
# 
# ggsave("figure_papier/figure_5.jpeg",
#        width = 12,
#        height = 12,
#        units = c("in")
# )
# 
# 
# 
# 
# write("\\end{document}", path_latex, append = TRUE)
# 
# 
# 
# ################################################################################
# # Supplementary materials
# path_latex_mat_spp <- "figure_papier/matsup.tex"
# 
# 
# df_pre_imput <- readRDS("data/genere/data_pre_impute_nb_vis_2.rds")
# 
# table_resume_pre_impute_rename <- df_pre_imput %>% rename_variables(var_instrum_name = var_instrumental_name)
# 
# 
# diag_var_table_pre_impute_latex <- table_resume_pre_impute_rename$complete_table %>%
#   mutate(diagnostic = !is.na(preffix_all)) %>%
#   select(Label, diagnostic)
# 
# 
# 
# 
# 
# 
# write(
#   "
# \\documentclass{article}
# \\usepackage[utf8]{inputenc}
# \\usepackage{float}
# \\usepackage{gensymb}
# \\usepackage{graphicx}
# \\usepackage{lscape}
# \\usepackage{longtable}
# \\usepackage{booktabs}
# \\usepackage{colortbl, xcolor}
# \\usepackage{hyperref}
# \\usepackage{subfiles}
# \\usepackage{setspace}
# \\usepackage{sectsty}
# \\usepackage{geometry}
# \\usepackage{indentfirst}
# \\usepackage{amsmath,lipsum}
# \\newcommand{\\mypm}{\\mathbin{\\smash{%
#   \\raisebox{0.35ex}{%
#     $\\underset{\\raisebox{0.2ex}{$\\smash -$}}{\\smash+}$%
#   }%
# }%
# }%
# }
# 
# \\usepackage{threeparttablex}
# \\usepackage{caption}                          
# \\DeclareCaptionLabelFormat{Sformat}{S#2 #1}    
# \\captionsetup[table]{labelformat=Sformat} 
# 
# 
# 
# \\usepackage[superscript,biblabel]{cite}
# \\sectionfont{\\clearpage}
# \\doublespacing 
# 
# \\title{Causal inference with multiple exposures: application of inverse-probability-of-treatment weighting to estimate the effect of daytime sleepiness in obstructive sleep apnea patients. \\\\
#   Supplementary Material}
# 
# \\author{}
# \\date{}
# 
# \\begin{document}
# 
# \\maketitle
# 
# \\subfile{authors}
# 
# ", path_latex_mat_spp,
# append = FALSE
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# write("
# \\clearpage
# \\section*{Supplementary material 1 : Imputation}
# Variables contained at least one missing which needed to be imputed, table \\ref{tab:Table_number_of_na} summarized number of missing values by adherence group for each of those variables.
# 
# We performed 10 imputed data sets using predictive mean matching as imputation method. After that we verified that the 10 imputed data sets were consistent with each other by looking at the distributions of the imputed variables in the different data sets.",
#       path_latex_mat_spp,
#       append = TRUE
# )
# 
# 
# 
# # Table des NA
# #  matériel supplémentaire
# nb_grp_obs <- max(select_var$data_frame$df_analyse$INS_obs_categ)
# df_latex_number_of_NA <-
#   df_pre_imput %>%
#   rbind(df_pre_imput %>%
#           mutate(INS_obs_categ = nb_grp_obs + 1)) %>%
#   group_by(INS_obs_categ) %>%
#   select_if(~ sum(is.na(.)) != 0) %>%
#   mutate(n = NA) %>% # tricks pour ajouter un compte des lignes
#   summarise_all(funs(sum(is.na(.)))) %>%
#   select(n, everything()) %>%
#   mutate_at(vars(-n, -INS_obs_categ), ~ paste0(formatC(., format = "f", big.mark = ",", digits = 0), " (", formatC(. * 100 / n, format = "f", big.mark = ",", digits = 1), "\\%)")) %>%
#   mutate(n = as.character(n)) %>%
#   t.df(pivot = "INS_obs_categ") %>%
#   rename_with(~ ifelse(as.numeric(. <= nb_grp_obs),
#                        paste0("Adherence group ", .),
#                        paste0("All adherence group")
#   ),
#   .cols = -key
#   ) %>%
#   left_join(table_rename$complete_table, by = c("key" = "var_name")) %>%
#   mutate(
#     Label = ifelse(is.na(real_name) & key == "n", "Number of patient", Label),
#     Label = ifelse(is.na(real_name) & key == "PPC_id_typeMasque", "Mask type", Label)
#   ) %>%
#   select(Label, `All adherence group`, contains("adherence group")) %>%
#   column_to_rownames(var = "Label")
# 
# df_latex_number_of_NA_order <- df_latex_number_of_NA %>%
#   rename_all(~ c(unlist(df_latex_number_of_NA["Number of patient", ], use.names = FALSE))) %>%
#   rownames_to_column() %>%
#   filter(rowname != "Number of patient") %>%
#   inner_join(diag_var_table_pre_impute_latex, by = c("rowname" = "Label")) %>%
#   mutate(rowname = R.utils::capitalize(str_remove(rowname, "^Diagnosis "))) %>%
#   arrange(desc(diagnostic))
# 
# 
# write(header_landscape, path_latex_mat_spp, append = TRUE)
# 
# df_latex_number_of_NA_order %>%
#   select(-diagnostic) %>%
#   mutate(
#     rowname = str_replace(rowname, "Epworth sleepiness scale", "ESS score"),
#     rowname = str_replace(rowname, "Gender", "Gender (male)"),
#     rowname = str_replace(rowname, "Apnea hypopnea index", " Apnea hypopnea index (event/h)")
#   ) %>%
#   kable(
#     format = "latex",
#     align = "c",
#     booktabs = TRUE,
#     escape = FALSE,
#     label = "Table_number_of_na",
#     linesep = "",
#     col.names = c("", colnames(df_latex_number_of_NA_order)[c(-1, -ncol(df_latex_number_of_NA_order))]),
#     caption = "Number of missing value by variables and group"
#   ) %>%
#   add_header_above(c("", "All groups", paste0(groupe_obs_table_latex, " (", seq_along(groupe_obs_table_latex), ")")),
#                    line = FALSE
#   ) %>%
#   pack_rows(
#     "Variables at diagnosis ", 1,
#     sum(df_latex_number_of_NA_order$diagnostic)
#   ) %>%
#   pack_rows(
#     "Variables at follow-up",
#     sum(df_latex_number_of_NA_order$diagnostic) + 1,
#     nrow(df_latex_number_of_NA_order)
#   ) %>%
#   kable_styling(latex_options = c(
#     "striped", "HOLD_position",
#     # "scale_down",
#     "repeat_header"
#   )) %>%
#   footnote(general = c("ESS : Epworth Sleepiness Scale; CPAP : Continuous Positive Airway Pressure"), general_title = "", threeparttable = TRUE) %>%
#   landscape() %>%
#   ecriture_fichier(path_latex_mat_spp)
# # write("\\clearpage", path_latex_mat_spp, append=TRUE)
# write(end_landscape, path_latex_mat_spp, append = TRUE)
# 
# write("\\clearpage", path_latex_mat_spp, append = TRUE)
# 
# 
# 
# 
# ##################################################################################
# 
# # troncature des poids
# summarry_boot_strap_pds_trunc %>%
#   arrondie_df(1) %>%
#   kable(
#     format = "latex",
#     align = "c",
#     booktabs = TRUE,
#     escape = FALSE,
#     label = "weight_truncations",
#     linesep = "",
#     caption = "Weight truncations"
#   ) %>%
#   kable_styling(latex_options = c("striped", "HOLD_position", "scale_down", "repeat_header")) %>%
#   footnote(general = c("Data are presented as mean (95% confidence interval) of bootstrap iterations"), general_title = "") %>%
#   ecriture_fichier(path_latex_mat_spp)
# write("\\clearpage", path_latex_mat_spp, append = TRUE)
# 
# #################################################################################
# ## Weighting models coefficients
# 
# weight_coef_boot_CI_rename <- table_rename$complete_table %>%
#   select(var_name, Label) %>%
#   right_join(weight_coef_boot_CI,
#              by = c("var_name" = "name")
#   ) %>%
#   select(-var_name) %>%
#   mutate(Label = str_replace(
#     Label, "Diagnosis apnea hypopnea index",
#     "Apnea Hypopnea Index at diagnosis"
#   )) %>%
#   mutate(
#     Label = str_replace(Label, "Epworth sleepiness scale", "ESS score"),
#     Label = str_replace(Label, "Gender", "Gender (male)"),
#     Label = str_replace(Label, "Apnea hypopnea index", " Apnea hypopnea index (event/h)")
#   )
# 
# weight_coef_boot_CI_rename %>%
#   mutate(coef = paste0(sprintf("%.03f", moy), "(", sprintf("%.03f", ci_inf), "; ", sprintf("%.03f", ci_sup), ")")) %>%
#   select(Label, rowname, coef) %>%
#   pivot_wider(names_from = rowname, names_sep = ".", values_from = coef) %>%
#   rename_all(~ c("Variables", paste0("Coefficients of ", groupe_obs_table_latex[-1], " group"))) %>%
#   kable(
#     format = "latex",
#     align = "c",
#     booktabs = TRUE,
#     escape = FALSE,
#     # longtable = TRUE,
#     label = "Wheight_coefficients",
#     linesep = "",
#     caption = "IPWRA weight model coefficient to assess the probability of being in an adherence group"
#   ) %>%
#   kable_styling(
#     font_size = 7,
#     latex_options = c(
#       "striped", "HOLD_position",
#       "repeat_header"
#     )
#   ) %>%
#   footnote(
#     general = c("Data are presented as mean (95% confidence interval) of bootstrap iterations", paste0(c("ESS : Epworth Sleepiness Scale", foot_note_var_presence), collapse = "; ")),
#     general_title = ""
#   ) %>%
#   ecriture_fichier(path_latex_mat_spp)
# write("\\clearpage", path_latex_mat_spp, append = TRUE)
# 
# 
# #################################################################################
# ## regression models coefficients
# ############### TEMP
# reg_coef_boot_CI_rename <- table_rename$complete_table %>%
#   select(var_name, Label) %>%
#   right_join(reg_coef_boot_CI,
#              by = c("var_name" = "name")
#   ) %>%
#   select(-var_name) %>%
#   mutate(Label = str_replace(
#     Label, "Diagnosis apnea hypopnea index",
#     "Apnea Hypopnea Index at diagnosis"
#   )) %>%
#   mutate(
#     Label = str_replace(Label, "Epworth sleepiness scale", "ESS score"),
#     Label = str_replace(Label, "Gender", "Gender (male)"),
#     Label = str_replace(Label, "Apnea hypopnea index", " Apnea hypopnea index (event/h)")
#   )
# 
# weight_coef_boot_CI_rename %>%
#   mutate(coef = paste0(sprintf("%.03f", moy), "(", sprintf("%.03f", ci_inf), "; ", sprintf("%.03f", ci_sup), ")")) %>%
#   select(Label, rowname, coef) %>%
#   pivot_wider(names_from = rowname, names_sep = ".", values_from = coef) %>%
#   rename_all(~ c("Variables", paste0("Coefficients of ", groupe_obs_table_latex[-1], " group"))) %>%
#   kable(
#     format = "latex",
#     align = "c",
#     booktabs = TRUE,
#     escape = FALSE,
#     # longtable = TRUE,
#     label = "Wheight_coefficients",
#     linesep = "",
#     caption = "IPWRA weight model coefficient to assess the probability of being in an adherence group"
#   ) %>%
#   kable_styling(
#     font_size = 7,
#     latex_options = c(
#       "striped", "HOLD_position",
#       "repeat_header"
#     )
#   ) %>%
#   footnote(
#     general = c("Data are presented as mean (95% confidence interval) of bootstrap iterations", paste0(c("ESS : Epworth Sleepiness Scale", foot_note_var_presence), collapse = "; ")),
#     general_title = ""
#   ) %>%
#   ecriture_fichier(path_latex_mat_spp)
# write("\\clearpage", path_latex_mat_spp, append = TRUE)
# 
# 
# 
# 
# ##################################################################################
# ## IPWRA models coefficients
# IPWRA_coef_boot_CI_rename <- table_rename$complete_table %>%
#   select(var_name, Label) %>%
#   right_join(IPWRA_coef_boot_CI %>% 
#                mutate(
#                  temp_num_name = str_extract(name, "\\d"),
#                  name = str_replace(name, "INS_obs_categ\\d$", "INS_obs_categ"),
#                  name = str_replace(name, "Vis_init_FDR_fumeur\\d$", "Vis_init_FDR_fumeur"),
#                  name = str_replace(name, "PPC_id_typeMasque", "PPC_id_typeMasque\\.")
#                ),
#              by = c("var_name" = "name")
#   ) %>%
#   mutate(
#     Label = str_remove(Label, "\\.\\d$"),
#     Label = ifelse(is.na(temp_num_name), Label, paste0(Label, " ", temp_num_name))
#   ) %>%
#   select(-var_name, -temp_num_name) %>%
#   mutate(Label = str_replace(
#     Label, "Diagnosis apnea hypopnea index",
#     "Apnea Hypopnea Index at diagnosis"
#   )) %>%
#   mutate(
#     Label = str_replace(Label, "Epworth sleepiness scale", "ESS score"),
#     Label = str_replace(Label, "Gender", "Gender (male)"),
#     Label = str_replace(Label, "Apnea hypopnea index", " Apnea hypopnea index (event/h)")
#   )
# 
# 
# IPWRA_coef_boot_CI_rename %>%
#   mutate(`Model coefficients` = paste0(sprintf("%.03f", moy), "(", sprintf("%.03f", ci_inf), "; ", sprintf("%.03f", ci_sup), ")")) %>%
#   select(Label, `Model coefficients`) %>%
#   kable(
#     format = "latex",
#     align = "c",
#     booktabs = TRUE,
#     escape = FALSE,
#     longtable = TRUE,
#     label = "IPWRA_coefficients",
#     linesep = "",
#     caption = "Multivariable weighted linear regression to assess the impact of CPAP adherence group on ESS under CPAP"
#   ) %>%
#   kable_styling(font_size = 6, latex_options = c(
#     "striped", "HOLD_position",
#     # "scale_down",
#     "repeat_header"
#   )) %>%
#   footnote(
#     general = c("Data are presented as mean (95% confidence interval) of bootstrap iterations", paste0(c("ESS : Epworth Sleepiness Scale", foot_note_var_presence), collapse = "; ")),
#     general_title = ""
#   ) %>%
#   ecriture_fichier(path_latex_mat_spp)
# write("\\clearpage", path_latex_mat_spp, append = TRUE)
# 
# 
# #################################################################################
# 
# 
# 
# 
# 
# 
# write("\\end{document}", path_latex_mat_spp, append = TRUE)
