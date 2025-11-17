

library("podacc")

data = amiodarone_tggates_29day_BMD

head(data)

list_bmd_values = data[,"BMD"]


results = runDistMethods(list_bmd_values)

# Examples how to plot results
plotDistMethodResults(list_bmd_values,results)

plotDistMethodResults(list_bmd_values,results,titleplot="DEMO Accumulation Plot Results",xlab_text = "BMD (mg/kg/day)")

plotDistMethodResults(list_bmd_values,results,titleplot="DEMO Accumulation Plot Results",xlab_text = "BMD (mg/kg/day)",GOBP=50)


