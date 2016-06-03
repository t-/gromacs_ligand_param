#!/bin/bash
fname=$1
sed -i 's/cag cag nhg hng  9   180.0       4.39320   2/cag cag nhg hng  9   180.0      15.16700   2; o_O modded - 4.39320 a488/g' $fname  #a488
sed -i 's/cag cwg cwg cag  9   180.0      10.66920   2/cag cwg cwg cag  9     0.0      10.66920   2; o_O modded - 180.0 a488/g' $fname  #a488
sed -i 's/cag cwg cwg cag  4  180.00       4.60240   2/cag cwg cwg cag  4    0.00       4.60240   2; o_O modded - 180.0 a488/g' $fname  #a488
sed -i 's/cag nhg ceg ceg  9   180.0       2.82420   2/cag nhg ceg ceg  9   180.0      27.82360   2; o_O modded - 2.82420 a647/g' $fname  #a647
sed -i 's/cag nhg c2g ceg  9   180.0       2.82420   2/cag nhg c2g ceg  9   180.0      27.82360   2; o_O modded - 2.82420 a647/g' $fname  #a647
sed -i 's/c2g ceg ceg cfg  9   180.0       4.18400   2/c2g ceg ceg cfg  9   180.0      27.82360   2; o_O modded - 4.18400 a647/g' $fname  #a647
sed -i 's/cfg ceg ceg nhg  9   180.0       4.18400   2/cfg ceg ceg nhg  9   180.0      27.82360   2; o_O modded - 4.18400 a647/g' $fname  #a647
sed -i 's/ceg cfg cfg ceg  9   180.0       4.18400   2/ceg cfg cfg ceg  9   180.0      27.82360   2; o_O modded - 4.18400 a647/g' $fname  #a647
sed -i 's/ceg ceg c2g nhg  9   180.0       4.18400   2/ceg ceg c2g nhg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy
sed -i 's/cfg cfg c2g nhg  9   180.0       4.18400   2/cfg cfg c2g nhg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy
sed -i 's/c3g c2g ceg ceg  9   180.0       4.18400   2/c3g c2g ceg ceg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy
sed -i 's/c3g c2g cfg cfg  9   180.0       4.18400   2/c3g c2g cfg cfg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy
sed -i 's/c3g ceg ceg cfg  9   180.0       4.18400   2/c3g ceg ceg cfg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy
sed -i 's/c3g cfg cfg ceg  9   180.0       4.18400   2/c3g cfg cfg ceg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy
sed -i 's/cfg ceg ceg nhg  9   180.0       4.18400   2/cfg ceg ceg nhg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy
sed -i 's/ceg cfg cfg nhg  9   180.0       4.18400   2/ceg cfg cfg nhg  9   180.0      27.82360   2; o_O modded cy/g' $fname #cy

echo -e '\E[35m made following modifications in '$(pwd)
cat $fname | grep modded
echo -e '\E[37m --stop no more modifications found'
