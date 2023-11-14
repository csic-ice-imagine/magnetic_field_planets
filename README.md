# magnetic_solar_planets


Repository used to reconstruct the planetary magnetic field with the usual Schmidt semi-normalised spherical harmonic coefficients for Earth, Jupiter, Saturn, Neptune and Uranus. All files from the directory data/ can be found (or constructed) with the open scientific data from the corresponding space missions.

To run the code in terminal you only need to run:

    python main.py

(or python3). Otherwise use anaconda, jupyter, etc.

To enable all the plotting option you need teh cartopy library, which ideally you will have enough by typing:

pip install cartopy

On some ubutnu version I had to fight a bit, and I got it by using:

sudo apt-get install libproj-dev proj-data proj-bin  
sudo apt-get install libgeos-dev  
sudo pip install cython  
sudo pip install cartopy

Free to use and change!
