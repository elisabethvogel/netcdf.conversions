# .onAttach <- function(libname, pkgname) { packageStartupMessage('Welcome to my package') } .onLoad <- function(libname,
# pkgname) { op <- options() op.devtools <- list( devtools.path = '~/R-dev', devtools.install.args = '', devtools.name =
# 'Elisabeth Vogel', devtools.desc.author = 'Elisabeth Vogel <elisabeth.vogel@climate-energy-college.org>', #
# devtools.desc.license = 'What license is it under?', devtools.desc.suggests = NULL, devtools.desc = list() ) toset <-
# !(names(op.devtools) %in% names(op)) if(any(toset)) options(op.devtools[toset]) invisible() }
