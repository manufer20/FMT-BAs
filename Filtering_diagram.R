# install.packages(c("DiagrammeR","DiagrammeRsvg","rsvg"))

library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

dot <- "
digraph consort_flow {
  graph [rankdir=TB, splines=ortho, pad=0.5];
  node  [shape=box, style=solid, fontname=\"Helvetica\", fontsize=10, margin=0.15];
  edge  [arrowhead=vee, color=\"#404040\"];

  /* Main vertical nodes */
  A [label=\"Initial Dataset\\n1070 Samples | 47 Donors\"];
  B [label=\"Post-Age Filter\\n777 Samples | 47 Donors\"];
  C [label=\"Post-Antibiotic Exclusion\\n583 Samples | 44 Donors\"];
  D [label=\"Post-Early Failure Removal\\n468 Samples | 42 Donors\"];
  E [label=\"Post-Frequency Filter\\n364 Samples | 17 Donors\"];
  F [label=\"Final Analytical Cohort\\n17 Donors\",
     shape=box,
     style=\"rounded,filled\",
     fillcolor=\"#E0F2F7\"] ;

  /* Side exclusions */
  X1 [label=\"Excluded for Age\\n(≤60 yrs)\\n293 Samples\"];
  X2 [label=\"Excluded for Antibiotic Use\\n194 Samples | 3 Donors\"];
  X3 [label=\"Excluded for Early Failure\\n115 Samples | 2 Donors\"];
  X4 [label=\"Excluded for Infrequent Donation\\n(<10 Samples)\\n104 Samples | 25 Donors\"];

  /* Invisible junction points */
  j1 [shape=point,width=0,height=0,label=\"\"];
  j2 [shape=point,width=0,height=0,label=\"\"];
  j3 [shape=point,width=0,height=0,label=\"\"];
  j4 [shape=point,width=0,height=0,label=\"\"];

  /* Flow arrows */
  A -> j1 [arrowhead=none];
  j1 -> B;
  j1 -> X1;

  B -> j2 [arrowhead=none];
  j2 -> C;
  j2 -> X2;

  C -> j3 [arrowhead=none];
  j3 -> D;
  j3 -> X3;

  D -> j4 [arrowhead=none];
  j4 -> E;
  j4 -> X4;

  E -> F;

  /* Align side boxes at same rank as junctions */
  {rank=same; j1; X1;}
  {rank=same; j2; X2;}
  {rank=same; j3; X3;}
  {rank=same; j4; X4;}
}
"

# Render the graph
gr <- grViz(dot)

# Export SVG and PNG
svg_text <- export_svg(gr)
cat(svg_text, file="consort_flow.svg")
rsvg_png(charToRaw(svg_text), file="consort_flow.png")

# View in RStudio / Notebook
print(gr)
browseURL("consort_flow.svg")
