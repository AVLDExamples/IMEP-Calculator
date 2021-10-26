# IMEP-Calculator
Example project to calculate IMEP from an AVL FIRE&amp;tm; full cylinder 4-stoke ICE engine simulation

Requiremens:
* AVL FIRE M R2021.1 or later (obviously)
* a finished simulation for a full 4 stroke cycle
* 2D output for the cylinder results:
  *  average pressure
  *  cylinder volume

Then, the IMEP/p_me  is calculated as$$ p_me = \int_0^720 p \,dV$$ and written to the results folder as a summary/KPI result
