#!/usr/bin/env Rscript
rm(list = ls(all = TRUE))

if (!require("pacman")) install.packages("pacman")
# list.of.packages <- c("tidyr","tidyverse", "magrittr", "reshape2", "dplyr", "RColorBrewer", "reshape2", "scales")
list.of.packages <- c("tidyverse", "magrittr")
pacman::p_load(list.of.packages, character.only = TRUE)

# ojo que açò són números de voxels
# X: 2480 um, 124 voxels de 20 um
# Y: 1680 um, 84 voxels de 20 um
# el centre són els voxels: 60, 61, 62, 63, 64, 65, però vull el endo en [65,68] -> pille [60-69]
ecm <- read.table("../config/setup_ecm.csv", header = F, sep = ",", stringsAsFactors = FALSE)
ecm2 <- ecm %>% filter(V1<69, 60<V1, V2<69, 60<V2) #està buit

blood <- read.table("../config/setup_endo_norm.csv", header = F, sep = ",", stringsAsFactors = FALSE)
blood2 <- blood %>% filter(V1<69, 60<V1, V2<69, 60<V2) %>% mutate(V1 = V1-60, V2 = V2-60) %>% slice(rep(1:n(), each = 8)) %>% mutate(V3 = rep(seq(1:8), 2))

# ecm3 <- ecm %>% filter((V1==0 & V2<9) | (V1==9 & V2<9))
# a1<-c(0,0,0)
# ecm3<-data.frame() %>% rbind(a1) %>% slice(rep(1:n(), each = 9)) %>% rename(V1 = X0, V2 = X0.1, V3 = X0.2)

N   <- 3
vec <- seq(1:10)-1
lst <- lapply(numeric(N), function(x) vec)
ecm2 <- expand.grid(lst) %>% as.matrix() %>% as.data.frame() %>% rename(V1=Var1, V2=Var2,V3=Var3)
ecm3 <- rbind(ecm2 %>% filter(V1==0 | V1==9),
              ecm2 %>% filter(V2==0 | V2==9),
              ecm2 %>% filter(V3==0 | V3==9),
              blood2
              )
ecm4a <- ecm3 %>% distinct(V1, V2, V3, .keep_all = TRUE)
ecm4 <- ecm3[!duplicated(ecm3[,1:3]),]

# a1<-ecm4 == ecm4a
# table(a1)["TRUE"] # 1464 = 488*3
# table(a1)["FALSE"] #NA

# volem les cèl.lules que estiguen en un quadrat de 200 um (10 voxels) com he pillat [60-69] i el centre es 64-65 això és (n-64)*20
cells <- read.table("../config/cells.csv", header = F, sep = ",", stringsAsFactors = FALSE)
# cells2 <- cells %>% filter(V1<100, -100<V1, V2<100, -100<V2)
cells2 <- cells %>% mutate(V1 = V1-(64*20/2), V2 = V2-(64*20/2)) %>% filter(-90<V1, V1<90, -90<V2, V2<90)
# range(cells2$V1)
# range(cells2$V2)
cells3 <- cells2 %>% slice(rep(1:n(), each = 8)) %>% mutate(V3 = rep(seq(1:8)*20-90, 30))
# range(cells3$V3)

write.table(blood2, file = "../config/endo_box.tsv",quote = F, row.names = F, col.names = F, sep="\t")
write.table(ecm4, file = "../config/ecm_box.tsv",quote = F, row.names = F, col.names = F, sep="\t")
write.table(cells3, file = "../config/cells_box.csv",quote = F, row.names = F, col.names = F, sep=",")
