
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/xin/Desktop/scalable/git/install/lib/

#gcc -o mytest timings2.c -I /home/xin/Desktop/scalable/p4est/install/include/ -I /usr/include/mpi/ -L /home/xin/Desktop/scalable/p4est/install/lib/ -lp4est -lmpi -lsc
#gcc -o 3dtest timings3.c -I /home/xin/Desktop/scalable/p4est/install/include/ -I /usr/include/mpi/ -L /home/xin/Desktop/scalable/p4est/install/lib/ -lp4est -lmpi -lsc

#gcc -c -o api.o api.c -I /home/xin/Desktop/scalable/p4est/install/include/ -I /usr/include/mpi/ 
#gcc -c -o cube.o cube.c  -I /home/xin/Desktop/scalable/p4est/install/include/ -I /usr/include/mpi/ 

gcc -c -o cube3.o cube3.c  -I /mnt/xin/install-mpi/include/ -I /usr/include/mpi/ -I /home/xin/Desktop/scalable/git/install/include
gcc -o cube3d cube3.o -L /home/xin/Desktop/scalable/git/install/lib/ -L /mnt/xin/install-mpi/lib/ -lp4est -lmpi -lsc  

gcc -c -o volume.o volume.c  -I /mnt/xin/install-mpi/include/ -I /usr/include/mpi/ -I /home/xin/Desktop/scalable/git/install/include
gcc -o volume volume.o -L /home/xin/Desktop/scalable/git/install/lib/ -L /mnt/xin/install-mpi/lib/ -lp4est -lmpi -lsc  


#gcc -o cube2d cube.o api.o -I /home/xin/Desktop/scalable/p4est/install/include/ -I /usr/include/mpi/ -L /home/xin/Desktop/scalable/p4est/install/lib/ -lp4est -lmpi -lsc
#echo $LD_LIBRARY_PATH
#sudo ldconfig
#./mytest -c rotcubes -l 1
