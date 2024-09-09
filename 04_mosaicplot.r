## mosaicplot

library(vcd)

struct <- structable(~ meta.cat+chic_risk_stratification, data = tumor)
mosaic(struct, , direction = "v", pop = FALSE,colorize = T, shade = TRUE,
       gp = gpar(fill = matrix(c("red","grey90" , "grey90","grey90" , "grey90", "green3"), 2, 3)))
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)