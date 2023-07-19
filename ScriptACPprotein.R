library(FactoMineR)
setwd("D:/FORMATION BIOIF ET BIOSTATISTIQUE/TRAVAIL DE GROUPE DIARRI/TD GROUPE GENETIQUE/ACP/ACP/Docs")
a=protein<- read.csv("protein.csv", header=TRUE,dec=".",row.names=1, sep=";",na.strings="")
colnames(protein)


library(ggplot2)
library(FactoMineR)
library(factoextra)
library(factoextra)
library(psych)
library(shiny)
library(FactoInvestigate)
library(Factoshiny)

library(missMDA)

library(FactoMineR)
library(factoextra)


res.pca=PCA(protein, scale.unit=TRUE, ncp=2, graph=T)
############# Exporter la figure dans un dossier sous format tiff
library(ggplot2)
library(missMDA)

eig.val <- get_eigenvalue(res.pca)
eig.val
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))

var <- get_pca_var(res.pca)
var
head(var$coord)
head(var$cos2)
head(var$contrib)
fviz_pca_var(res.pca, col.var = "black")

Qualité de représentation

La qualité de représentation des variables sur la carte de l’ACP s’appelle cos2 (cosinus carré) . Vous pouvez accéder au cos2 comme suit:
  
  head(var$cos2, 4)

#Vous pouvez visualiser le cos2 des variables sur toutes les dimensions en utilisant le package corrplot:
  
  library("corrplot")
## corrplot 0.92 loaded
corrplot(var$cos2, is.corr=FALSE)


#####Il est également possible de créer un bar plot du cosinus carré des variables en utilisant la fonction fviz_cos2() [dans factoextra]:
  
  # Cos2 total des variables sur Dim.1 et Dim.2
  fviz_cos2(res.pca, choice = "var", axes = 1:2)

  # Colorer en fonction du cos2: qualité de représentation
  fviz_pca_var(res.pca, col.var = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE # évite le chevauchement de texte
  )
  head(var$contrib, 4)
  library("corrplot")
  corrplot(var$contrib, is.corr=FALSE)   
  # Contributions des variables ? PC1
  fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
  
  
  # Contributions des variables ? PC2
  fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
  
  
  La contribution totale ? PC1 et PC2 est obtenue avec le code R suivant:
    
    fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 10)
  
  library(missMDA)
  #Graphique des individus Résultats
  
 # Les résultats, pour les individus, peuvent être extraits à l’aide de la fonction get_pca_ind() [package factoextra]. Comme get_pca_var(), la fonction get_pca_ind() retourne une liste de matrices contenant tous les résultats pour les individus (coordonnées, corrélation entre individus et axes, cosinus-carré et contributions)
  
  ind <- get_pca_ind(res.pca)
  ind
  # Coordonnées des individus
  head(ind$coord) 
  # Qualit? des individus
  head(ind$cos2)
  # Contributions des individus
  head(ind$contrib)  
  fviz_pca_ind (res.pca)  
  
  library(missMDA)
  #Comme les variables, il est également possible de colorer les individus en fonction de leurs valeurs de cos2:
    
  fviz_pca_ind (res.pca, col.ind = "cos2",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE # évite le chevauchement de texte
  )

  Notez que les individus qui sont similaires sont regroupés sur le graphique.
  
  fviz_pca_ind (res.pca, pointsize = "cos2",
                pointshape = 21, fill = "#E7B800",
                repel = TRUE # ?vite le chevauchement de texte
  )                  

  fviz_pca_ind(res.pca, col.ind = "cos2", pointsize = "cos2",
               gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
               repel = TRUE 
  ) 
  
  # Biplot of individuals and variables
  fviz_pca_biplot(res.pca, repel = TRUE)
  
  fviz_cos2(res.pca, choice = "ind")

  # Contribution totale sur PC1 et PC2
  fviz_contrib(res.pca, choice = "ind", axes = 1:2)  
 
  
  Cas de variables continues
  
  1 Faire une ACP
  2 Appliquer la classification hiérarchique sur le résultat de l’ACP
  library(FactoMineR)
  # 1. ACP 
  res.pca=PCA(protein, scale.unit=TRUE, ncp=3, graph=F)
  # 2. HCPC
  res.hcpc <- HCPC(res.pca, graph = FALSE)
  3 Visualiser le dendrogramme généré par la classification. Fonction R: fviz_dend() [factoextra]:
    fviz_dend(res.hcpc, 
              cex = 0.7,                     # Taille du text
              palette = "jco",               # Palette de couleur ?ggpubr::ggpar
              rect = TRUE, rect_fill = TRUE, # Rectangle autour des groupes
              rect_border = "jco",           # Couleur du rectangle
              labels_track_height = 0.8      # Augment l'espace pour le texte
    )
  

  
  Le dendrogramme suggére une solution à 3 groupes.
  
  Visualiser les individus et colorer par groupes. Fonction R: fviz_cluster() [factoextra].
  
  fviz_cluster(res.hcpc,
               repel = TRUE,            # Evite le chevauchement des textes
               show.clust.cent = TRUE, # Montre le centre des clusters
               palette = "jco",         # Palette de couleurs, voir ?ggpubr::ggpar
               ggtheme = theme_minimal(),
               main = "Factor map"
  )
  
  
  
  Contenu du résultat de la fonction HCPC():
    
    data.clust: Données d’origine avec une colonne supplémentaire appelée clust contenant les groupes.
  desc.var: les variables décrivant les groupes
  desc.ind: les individus les plus typiques de chaque groupes
  desc.axes: les axes décrivant les groupes
  Données d’origine avec la colonne class:
    
    head(res.hcpc$data.clust, 10)
 
  Variables quantitatives décrivant le plus chaque cluster:
    
    res.hcpc$desc.var$quanti
    
    r=res.hcpc$desc.var$quanti
  
  Les résultats ci-dessus indiquent que les individus dans le groupes 1 ont des coordonnées élevées sur l’axes 1 et sur l’axe 2. Les individus du groupe 2 ont des coordonnées élevées sur le deuxiéme axe. Les individus appartenant au troisiéme groupe ont des coordonnées élevées sur l’axe 1.
  
  Individus représentatifs de chaque groupe:
    
    res.hcpc$desc.ind$para
  
##################################################
----
    