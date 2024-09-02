# Loading packages --------------------------------------------------------
my_dir = getwd()
setwd(my_dir)
### Packages / library ###

library(readxl)
library(ggplot2)
library(dplyr)
library(grDevices)


# Data loading body positions------------------------------------------------------------

ECP_body_position <- read_excel("Data_ECP_final_body_position.xlsx")

# Data formatting body positions ---------------------------------------------------------

Data_to_analyse = data.frame(Body_pos=character(),ECP=double())

row_number = nrow(ECP_body_position)

for (cpt in 1:row_number){
  if (ECP_body_position$Temps_face_face[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("face-to-face",ECP_body_position$ECP_face_face[cpt])
  } 
  if (ECP_body_position$Temps_face_dos[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("face-to-back",ECP_body_position$ECP_face_dos[cpt])
  } 
  if (ECP_body_position$Temps_face_flanc[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("face-to-flank",ECP_body_position$ECP_face_flanc[cpt])
  } 
  if (ECP_body_position$Temps_face_bd[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("face-to-dp",ECP_body_position$ECP_face_bd[cpt])
  } 
}

# Graphic body positions -----------------------------------------------------------------

ggplot(Data_to_analyse,aes(as.factor(Body_pos),ECP)) +
  geom_boxplot(fill=c("grey", "grey","green","grey"))+
  theme_classic()+ theme(legend.position = "none")

wilcox.test(Data_to_analyse$ECP[Data_to_analyse$Body_pos=='face-to-face'],Data_to_analyse$ECP[Data_to_analyse$Body_pos=='face-to-flank'])

# Data loading body parts------------------------------------------------------------

ECP_body_parts <- read_excel("C:/Thèse/Rédaction/Article Legrand et al. 2024/Pour soumission/Raw data/ECP/Data_ECP_final_body_parts.xlsx")

Data_to_analyse = data.frame(Body_part=character(),ECP=double())

row_number = nrow(ECP_body_parts)

for (cpt in 1:row_number){
  if (ECP_body_parts$Temps_face[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("face",ECP_body_parts$ECP_face[cpt])
  } 
  if (ECP_body_parts$Temps_back_head[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("backhead",ECP_body_parts$ECP_back_head[cpt])
  } 
  if (ECP_body_parts$Temps_neck[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("neck",ECP_body_parts$ECP_neck[cpt])
  } 
  if (ECP_body_parts$Temps_ventre[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("ventre",ECP_body_parts$ECP_ventre[cpt])
  } 
  if (ECP_body_parts$Temps_chest[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("chest",ECP_body_parts$ECP_chest[cpt])
  } 
  if (ECP_body_parts$Temps_vp[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("vp",ECP_body_parts$ECP_vp[cpt])
  } 
  
  if (ECP_body_parts$Temps_flank[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("flank",ECP_body_parts$ECP_flank[cpt])
  } 
  
  if (ECP_body_parts$Temps_dp[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("dp",ECP_body_parts$ECP_dp[cpt])
  } 
  
  if (ECP_body_parts$Temps_back[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("back",ECP_body_parts$ECP_back[cpt])
  } 
  
  if (ECP_body_parts$Temps_members[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("members",ECP_body_parts$ECP_memebers[cpt])
  } 
  
  if (ECP_body_parts$Temps_Armpit[cpt] > 0){
    nrow_temp = nrow(Data_to_analyse)
    Data_to_analyse[nrow_temp+1,] = list("armpit",ECP_body_parts$ECP_armpit[cpt])
  } 
}

# Graphic body parts -----------------------------------------------------------------
#fill=c("grey", "grey","green","grey")
ggplot(Data_to_analyse,aes(as.factor(Body_part),ECP)) +
  geom_boxplot(fill=c("grey", "grey","grey","green","grey","green","grey","grey","green","grey","grey"))+
  theme_classic()+ theme(legend.position = "none")

