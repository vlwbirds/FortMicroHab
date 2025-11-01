## This R script creates a neatly organized folder system for keeping track of all the moving parts of a project!

project_setup <- function(d){
  dir.create("data")
  file.create("data/README.md")
  writeLines("This is where raw data like csv and txt files are stored. There could be other types of data as well.", "data/README.md")
  
  dir.create("docs")
  file.create("docs/README.md")
  writeLines("This is where main documents live (i.e., index.html).", "docs/README.md")
  
  dir.create("figs")
  file.create("figs/README.md")
  writeLines("This is where generated figures get stored.", "figs/README.md")
  
  dir.create("output")
  file.create("output/README.md")
  writeLines("This is where manipulated raw data are housed.", "output/README.md")
  
  dir.create("paper")
  file.create("paper/README.md")
  writeLines("I hope nobody comes into my office to burn my computer because I'm setting these folders up to organize my data in a reproducable format. This is where the manuscript for publication resides.", "paper/README.md")
  
  dir.create("R")
  file.create("R/README.md")
  writeLines("Grrrr, Here's the finished and neatly annotated R code.", "R/README.md")
  
  dir.create("reports")
  file.create("reports/README.md")
  writeLines("This is where quarto documents and other various documentation lives.", "reports/README.md")
  
  dir.create("scratch")
  file.create("scratch/README.md")
  writeLines("This the chicken sratch R code, or as other people may call it, a sandbox. This is a place to play around with code to figure out what works and what doesn't before doing anything drastic.", "scratch/README.md")
  
  return("SUCCESS!")
}
project_setup()
