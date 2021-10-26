# IMEP-Calculator
Example project to calculate IMEP from an AVL FIRE&amp;trade; full cylinder 4-stoke ICE engine simulation

Requiremens:
* AVL FIRE M R2021.1 or later (obviously)
* a finished simulation for a full 4 stroke cycle
* 2D output for the cylinder results:
  *  average pressure
  *  cylinder volume

Then, the IMEP/p_me  is calculated as  <img src="https://render.githubusercontent.com/render/math?math=p_{me}=\int_0^720p\,dV"/> and written to the results folder as a summary/KPI result
