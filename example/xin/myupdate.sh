
cd /mnt/xin/p4est/
make

echo "...finish compile"

rm -f /mnt/xin/myp4est.tar.gz
if [ ! -f /mnt/xin/myp4est.tar.gz ]; then
    #echo "File not found! compressing..."
    cd /mnt/xin/
    tar zcf myp4est.tar.gz ./p4est
    tar zcf myinstall.tar.gz ./install-p4est
    #tar zcf /mnt/xin/spark.tar.gz ./spark-1.5.0
fi

echo "...finish compressing and start loop"
declare -a arr=("kid112" "kid113" "kid115")

#scp /mnt/xin/myp4est.tar.gz kid112:/mnt/xin/ 
#declare -a arr=("maquis11")

## now loop through the above array
for i in "${arr[@]}"
do
   echo "$i"
   ssh $i rm -rf /mnt/xin/p4est
   ssh $i rm -rf /mnt/xin/install-p4est
   scp /mnt/xin/myp4est.tar.gz $i:/mnt/xin/
   scp /mnt/xin/myinstall.tar.gz $i:/mnt/xin/
   #ssh $i tar zxf /mnt/xin/myp4est.tar.gz -C /mnt/xin/ 
   ssh $i "cd /mnt/xin; tar zxf myp4est.tar.gz"
   ssh $i "cd /mnt/xin; tar zxf myinstall.tar.gz"
done
