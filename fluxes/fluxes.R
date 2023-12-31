library(tidyverse)
library(fluxweb)

att = readRDS("fluxes/attributes_imputed.RData")

mat = as.matrix(read.csv("fluxes/metafoodweb.csv", row.names = 1))

thousand = vector(mode = "list",length=1000)
set.seed(404)
ind = rep(1:100, 10)

for (k in 1:length(ind)) { 
  
  # a list for feeding preference matrices
  mat.prefs = vector(mode = "list",
                     length=length(att))
  
  web = vector(mode = "list", 
               length=length(att))
  
  for (i in 1:length(att)) {
    web[[i]] = mat[att[[i]][[ind[k]]]$taxon,
                   att[[i]][[ind[k]]]$taxon]
  }
  
  for (i in 1:length(att)) {
    ####################   Omnivores' Balanced Diet Plan   #########################
    # add biomass values in the matrix to 'manually' define the preferences
    # first create a matrix with species biomasses
    mat.bioms = replicate(length(att[[i]][[ind[k]]]$Biomass.mg), att[[i]][[ind[k]]]$Biomass.mg)
    # mat.prefs contains preference of predators based on their prey biomasses
    mat.prefs[[i]] = web[[i]] * mat.bioms
    
    basals = which(att[[i]][[ind[k]]]$taxon %in% c("roots","detritus","bacteria","fungi"))
    animals = which(!(att[[i]][[ind[k]]]$taxon %in% c("roots","detritus","bacteria","fungi")))
    #omnivores that feed on basals and animals
    omnivores = which(colSums(mat.prefs[[i]][basals,])>0 &
                        colSums(mat.prefs[[i]][animals,])>0)
    # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
    mat.prefs[[i]][animals, omnivores] = mat.prefs[[i]][animals, omnivores, drop=FALSE] %*%
      diag(1/colSums(as.matrix(mat.prefs[[i]][animals, omnivores, drop=FALSE])),
           length(omnivores),length(omnivores)) #diag(4)!=diag(4,1,1) important if single omnivore
    
    # now we additionally make this sum to the complement of whatever else they eat
    std = colSums(web[[i]][animals, 
                           omnivores])
    mat.prefs[[i]][animals, 
                   omnivores] = mat.prefs[[i]][animals, 
                                               omnivores]*rep(std, 
                                                              each=nrow(mat.prefs[[i]][animals, 
                                                                                       omnivores]))
    
    mat.prefs[[i]] = vegan::decostand(mat.prefs[[i]], "total", MARGIN = 2)
  }
  
  # a list for flux matrices
  fluxes = vector(mode = "list",length=length(att))
  
  allmetrics = data.frame(Plot = character(length(att)),
                          Treatment = character(length(att)),
                          tot.flux = numeric(length(att)),   # total energy flux
                          pred.flux = numeric(length(att)),  # predation flux
                          herb.flux = numeric(length(att)),  # herbivory flux
                          detr.flux = numeric(length(att)),
                          secon.decomp.flux = numeric(length(att)),
                          down = numeric(length(att)),   # from herbivores per unit herbivore biomass
                          up = numeric(length(att)),     # to herbivores per unit herbivore biomass
                          herb.press = numeric(length(att))) # to herbivores per unit plant biomass
  
  
  for (i in 1:length(att)) {
    
    ################################# Uncertainty ##################################
    # Here we take each consumer in the foodweb and replace its fixed preferences 
    # with a random sample from a dirichlet distribution whose component probabilities 
    # are given by the vector of the original preferences. The vector is multiplied 
    # by a scalar that modifies the shape of the distribution (larger = less uncertainty)
    # Across several iterations our expectations regarding what consumers feed on 
    # are met, on average. But in each iteration consumer preferences deviate somewhat
    # from those expected based on intrinsic preference and/or relative availability.
    for (j in 1: nrow(mat.prefs[[i]])) { 
      mat.prefs[[i]][,j] = LaplacesDemon::rdirichlet(1, mat.prefs[[i]][,j]*100)
    }
    mat.prefs[[i]][is.nan(mat.prefs[[i]])] = 0 #removes NaNs from basal node "preferences"
    ################################################################################
    
    # down-weighing cannibalism
    diag(mat.prefs[[i]]) = diag(mat.prefs[[i]])*0.01
    
    fluxes[[i]] <- fluxing(mat.prefs[[i]],
                           att[[i]][[ind[k]]]$Biomass.mg, 
                           att[[i]][[ind[k]]]$Pop.met.rate.J_h,
                           att[[i]][[ind[k]]]$efficiency,
                           bioms.prefs = F,
                           bioms.losses = F,
                           ef.level = "prey")
    
    animals = which(!(att[[i]][[ind[k]]]$taxon %in% c("roots","detritus","bacteria","fungi")))
    plants = which(att[[i]][[ind[k]]]$taxon == "roots")
    detritus = which(att[[i]][[ind[k]]]$taxon == "detritus")
    microbs = which(att[[i]][[ind[k]]]$taxon %in% c("bacteria","fungi"))
    
    herbivores = which(colSums(fluxes[[i]][c(animals,
                                             detritus,
                                             microbs),,drop = FALSE]) == 0 &
                         colSums(fluxes[[i]][plants,,drop = FALSE])>0)
    
    plant.consumers = which(colSums(fluxes[[i]][plants,,drop = FALSE])>0)
    
    # the proportion of energy uptake in plant consumers that comes from plants
    # influx of energy from plants / total influx of energy
    # (1 for herbivores, a fraction for omnivores)
    prop = (fluxes[[i]][plants,
                        plant.consumers]*0.5446309)/
      colSums(fluxes[[i]][,plant.consumers]*att[[i]][[ind[k]]]$efficiency)
    
    predators = which(colSums(fluxes[[i]][c(plants,
                                            detritus,
                                            microbs),,drop = FALSE]) == 0 &
                        colSums(fluxes[[i]][animals,]) > 0)
    
    allmetrics[i,]$Plot = unique(att[[i]][[ind[k]]]$Plot)
    allmetrics[i,]$Treatment = unique(att[[i]][[ind[k]]]$Treatment)
    allmetrics[i,]$tot.flux = sum(fluxes[[i]])                         # total energy flux              
    allmetrics[i,]$pred.flux = sum(fluxes[[i]][animals, ])             # predation flux
    allmetrics[i,]$herb.flux = sum(fluxes[[i]][plants, ])              # herbivory flux
    allmetrics[i,]$detr.flux = sum(fluxes[[i]][detritus, ])            # detritivory flux
    allmetrics[i,]$secon.decomp.flux = sum(fluxes[[i]][microbs, ])     # secondary decomposers flux
    
    allmetrics[i,]$down = sum(fluxes[[i]][plant.consumers, ]*prop)    # outflux from plant consumers***
    allmetrics[i,]$up = sum(fluxes[[i]][plants, ])*0.5446309          # influx from plants to plant consumers
    
    
  }
  
  thousand[[k]] = allmetrics
  
  ############################## Show loop progress ##############################
  cat('\014')
  cat(paste0(round((k/1000) * 100), '% completed'))
  Sys.sleep(.05)
  if (k == 1000) cat(': Done')
  ################################################################################ 
}


thousand = do.call(rbind, thousand)

write.csv(thousand, "fluxes/fluxes1000_dir100.csv",
          row.names = F)
