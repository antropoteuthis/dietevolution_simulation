# Analayzing the Organismal and Evolutionary Determinants of Diet

Goals:

	-Simulate the determination of diet in a clade of predators with evolving functional traits
	-Compare the effects of different numbers of organismal determinants of diet
	-Compare the effects of different evolutionary histories of the organismal determinants of diet
	-Explore the information (general and evolutionary) contained in the unexplained variance of the diet 
  (not explained by measured determinants or measured prey-availability)

Simulation:

The pipeline uses real positions in Z (depth), relative abundances, and patchiness (degree of spatial aggregation, when available from transects) from the VARS data, and generates random X and Y components over thousands of iterations. In each iteration, my pipeline will record the relative frequencies of 2 point objects being within threshold (i.e. 1m) distance for each siphonophore-prey taxon couple. The output would be a relative encounter rate matrix for each siphonophore species (rows) and prey taxon (columns) in a scale from 0 to 1. Adjustments on the threshold distance can be made for different siphonophore and prey species, accounting for siphonophore colony length, tentacle length, prey motility, and prey size.

	-Simulated tree
	-Different numbers of “measured” and “unmeasured” traits simulated on the tree.
	-Selectivity indices determined for each prey type and species as a product of the deviations of different subsets of those evolving traits. 
	-Relative availabilities of prey items per species generated at random in 4 habitat-zones. Assingment of species to these habitat-zones depends on an evolving component and a random component.
	-Diet calculated as the product of preyfield and selectivity.

Analysis:

	-Mantel testing (regular and partial) for congruence signal in the diet.
	-Partitioning the unexplained variance, how much of it is phylogenetically structured, how much of it is not (random errors, measurement errors).
	-Assessing the number and the distributions of unmeasured traits.
	-Assess the evolutionary generative model behind the unmeasured traits.
