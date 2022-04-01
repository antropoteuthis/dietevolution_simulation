# Analayzing the Organismal and Evolutionary Determinants of Diet

I simulated the evolution of characters on the tree, including a habitat vector with a stochastic component. I generated random prey availabilities for each species depending on their habitat. In addition, I used linear models that would combine and transform evolving characters into specific selectivity loadings (mechanistic rules). Diet was calculated as a matrix multiplication between selectivity and availability. A variable number of observable and unobservable traits were simulated under different models of evolution (BM, OU, EB).  The resulting data was analyzed using variance partitioning tools, distance based methods, and phylogenetic methods. The results show the ability of the proposed methods to recover the relative contribution of the different generative processes in the diet-habitat-morphology-phylogeny data, and to identify the evolutionary models underlying the characters involved. Moreover, I was able to use the phylogenetic structure in the residuals to separate the relative contribution of unmeasured key traits from other sources of error.

Goals:

	-Simulate the determination of diet in a clade of predators with evolving functional traits
	-Compare the effects of different numbers of organismal determinants of diet
	-Compare the effects of different evolutionary histories of the organismal determinants of diet
	-Explore the information (general and evolutionary) contained in the unexplained variance of the diet 
  (not explained by measured determinants or measured prey-availability)

Simulation:

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
