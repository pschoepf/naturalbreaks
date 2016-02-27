# naturalbreaks
Ports of fisher/jenks natural breaks algorithm with O(k×n×log(n)) complexity, invented and described by Maarten Hilferink:
http://wiki.objectvision.nl/index.php/Fisher%27s_Natural_Breaks_Classification  

The code provided is a (nearly 1:1, as there is not much left for optimization) port of the corresponding C code that was
originally developed by Maarten Hilferink, © Object Vision BV and can be found here:
http://wiki.objectvision.nl/index.php/CalcNaturalBreaksCode

My contribution is:
- a basic port to Java
- (coming soon) a basic port to javascript so that the algorithm can be directly used in the browser to compute natural looking color scales.


License is provided on the same GNU GPL v3.0 license as the original contribution from Marteen.
