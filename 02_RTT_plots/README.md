# Distribution of the sea spider fossil record through time

The in-house R script [`Create_RTT_plot.R`](scripts/Create_RTT_plot.R) is used to generate the range through time plot that shows the distribution of the sea spider fossil record. The function used to create the plot modifies an original function within the R package `paleoverse` so that the range plot can be ordered by groups above species. The in-house function can be found in [`Functions_RTTplot.R](scripts/Functions_RTTplot.R).

If you run the [`Create_RTT_plot.R` script](scripts/Create_RTT_plot.R), you will obtain the same plot that you can find inside the [`plots` directory](plots/Pycnogonid_fossils.pdf). The data used were retrieved from the [The Paleobiology Database](http://paleobiodb.org/data1.2/occs/list.csv?datainfo&rowcount&base_name=Pycnogonida&ident=all&pgm=gplates,scotese,seton&show=full,ref), accessed on 2024-07-18 at 19:47:21 GM.
