# Modified SIMLANDER Code

This folder contains a modified version of the `allocator.R` script from the SIMLANDER v2.1.1 model.

To facilitate the changes inside the code search for the varibale **model_accessibility** (it allows you to change the weight of road and transit network) and the word **zoning** (it will briefly explain how it works). 
## 📌 Changes Made:
- Integrated **transit accessibility** as an additional factor in the urban allocation process.
- Adjusted weight calculation to prioritize pixels closer to public transport.
- Aligned land-use inputs and output formats to match the TOD scenario framework used in the thesis.

## 🔁 Based on:
Original SIMLANDER code by Richard Hewitt  
Version: [v2.1.1](https://simlander.wordpress.com/2025/02/20/simlander-v2-1-1-released/)  
License: Free and Open Source (as described by author; see thesis references)

## 📚 Reference
For detailed context, see Chapter 4 of the thesis:
> Dimas, P. (2025). *Modeling Transit-Oriented Development: A Land-Use Approach for Urban Decarbonization*. Master’s Thesis, University of Twente.
