![](./aux/TKSA-MC.jpg)
## Tanford-Kirkwood Surface Accessibility - Monte Carlo (TKSA-MC)

The TKSAMC is a web server which calculates protein charge-charge interactions via the Tanford-Kirkwood Surface Accessibility model with the Monte Carlo method for sampling different protein protonation states. 

The web server is freely available at UNESP (São Paulo State University - DF/IBILCE): http://tksamc.df.ibilce.unesp.br.

Example:
```bash
python tksamc.py -f sample_pdb_1ubq.pdb -ph 7.0 -T 300.0
```

>"TKSA-MC: A Web Server for rational mutation through the optimization of protein charge interactions -  
Vinícius G Contessoto, Vinícius M de Oliveira, Bruno R Fernandes, Gabriel GSlade, Vitor B. P. Leite, bioRxiv 221556; doi: https://doi.org/10.1101/221556"
