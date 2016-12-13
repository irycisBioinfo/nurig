library(dplyr)
library(tidyr)
library(XML)

#lista = read.table("lista.txt")
#colnames(lista) = c("file")
#ref = "b2.ref"
#outFile = "out.xml"



j = 1
for (i in lista$file)
{
  system(paste("nucmer ", ref, i, "> /dev/null", collapse = ""))
  system("delta-filter -q out.delta > out.f.delta ")
  system("show-coords -l -L 200 -q -B out.f.delta > out.coords")
  tmpTable = read.table("out.coords", sep = "\t")
  colnames(tmpTable) = c(
    "qID",
    "date",
    "lenq",
    "alitype",
    "referencefile",
    "rID",
    "startquery",
    "endquery",
    "startreference",
    "endreference",
    "peridentity",
    "persimilarity",
    "lenAlig",
    "c16",
    "c17",
    "c18",
    "c19",
    "strand",
    "lenref",
    "c21",
    "c22"
  )
  tmpTable$slot = j
  tmpTable$file = i
  
  if (j > 1)
  {
    finalTable = bind_rows(finalTable, tmpTable)
  } else{
    finalTable = tmpTable
  }
  j = j + 1
}

finalTable = finalTable %>% group_by(rID) %>% summarise(ln = first(lenref)) %>% mutate(prueba = cumsum(ln)) %>% mutate(position = prueba -
                                                                                                                         ln) %>% select(rID, position) %>% full_join(finalTable, .)
finalTable = finalTable %>% group_by(slot) %>% mutate(orden = sum(lenAlig *
                                                                    peridentity / 100)) %>% arrange(desc(orden))



tmp = finalTable %>%  split(.$slot)
size = finalTable %>% group_by(rID) %>% summarise(tl = first(lenref)) %>% summarise(S = sum(tl))
xml = xmlTree()
xml$addTag(
  "cgview",
  attrs = c(backboneRadius = "160", sequenceLength = as.integer(size[1])),
  close = FALSE
)



j = 1

if(palette == 'rainbow')
{
  colores = rainbow(max(finalTable$slot))
}else if (palette == 'topo'){
  colores = topo.colores(max(finalTable$slot))
}else if (palette == 'terrain'){
  colores = terrain.colors(max(finalTable$slot))
}

for (tab in 1:length(tmp))
{
  xml$addTag("featureSlot",
             attrs = c(strand = "direct"),
             close = FALSE)
  clr = paste("rgb(", paste(col2rgb(colores[j]), collapse = ","), ")", collapse = "")
  j = j + 1
  xml$addTag("feature",
             attrs = c(color = clr, decoration = "arc"),
             close = FALSE)
  for (r in 1:nrow(tmp[[tab]]))
  {
    #xml$addTag("feature", close = FALSE)
    xml$addTag(
      "featureRange",
      attrs = c(
        start = tmp[[tab]]$startreference[r] + tmp[[tab]]$position[r],
        stop = tmp[[tab]]$endreference[r] + tmp[[tab]]$position[r],
        opacity = tmp[[tab]]$peridentity[r] / 100
      ),
      close = FALSE
    )
    xml$closeTag()
    
  }
  xml$closeTag()
  xml$closeTag()
  
  
}
xml$closeTag()

xml$addTag("legend",
           attrs = c(position = "upper-right"),
           close = FALSE)
j = 1

tmp2 = finalTable %>% group_by(file, orden) %>% summarise(slot = first(slot)) %>% arrange(desc(orden))

for (i in tmp2$file)
{
  clr = paste("rgb(", paste(col2rgb(colores[j]), collapse = ","), ")", collapse = "")
  xml$addTag(
    "legendItem",
    attrs = c(
      text = i,
      drawSwatch = "true",
      swatchColor = clr
    ),
    close = TRUE
  )
  j = j + 1
}
xml$closeTag()

##### Dibujar features (contigs)

contigs = finalTable %>% group_by(rID) %>% summarise(position = first(position) +
                                                       1)
if (nrow(contigs) > 1)
{
  xml$addTag("featureSlot",
             attrs = c(strand = "direct"),
             close = FALSE)
  xml$addTag("feature",
             attrs = c(decoration = "clockwise-arrow"),
             close = FALSE)
  cl = rep(c("red", "blue"), times = (nrow(contigs) / 2) + 1)
  for (i in 1:(nrow(contigs) - 1))
  {
    xml$addTag(
      "featureRange",
      attrs = c(
        start = contigs$position[i],
        stop = contigs$position[i + 1],
        mouseover = contigs$rID[i],
        color = cl[i],
        showLabel = "FALSE"
      ),
      close = TRUE
    )
  }
  xml$addTag(
    "featureRange",
    attrs = c(
      start = contigs$position[nrow(contigs)]
      ,
      stop = as.integer(size[1]),
      mouseover = contigs$rID[nrow(contigs)],
      color = cl[i + 1],
      showLabel = "FALSE"
    ),
    close = TRUE
  )
  xml$closeTag()
  xml$closeTag()
}


saveXML(xml$value(), file = "out.xml")

