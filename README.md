# Gene Annotator
The gene-annotator is written in python3 that takes a list of genes and returns information from select databases.

## Supported Websites
* [Microbesonline](http://www.microbesonline.org)
* [NCBI Conserved Domain](https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi)


## Set up
* Check for access to websites (see IP Address Issues)
* Clone git repository: https://github.com/biofilmDB/gene-annotator.git
    
## Run
* Create the conda environment
    * conda env create -f environment.yml 
* run: 
        python driver.py your-gene-list-name.txt desired-output-file-name.csv

        
## Issues
* Microbesonline.org
    * Whitelisting IP to connect
        * To test connection run: 
                        mysql -h pub.microbesonline.org -u guest -pguest genomics
        * [Instructions to whitelist your IP](http://www.microbesonline.org/programmers.html#connectsql)
    * ConnectionResetError: [Errno 104] Connection reset by peer

## Contact
If you have any suggestions please contact [britney.gibbs@student.montana.edu](mailto:britney.gibbs@student.montana.edu)
