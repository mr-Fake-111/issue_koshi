import java.io.File
import java.io.FileWriter
import kotlin.math.max
import kotlin.math.pow
import kotlin.math.sqrt

class koshiSolver {

    companion object {

        fun kutta2ConstCheck(syst: (DoubleArray) -> DoubleArray, startVals: DoubleArray, end:Double, e: Double, h_ch: Double = 0.0): MutableList<DoubleArray> {

            val delta = (1.0/max(startVals[0], end)).pow(3) + norm2(syst(startVals)).pow(3) //выбор оптимального шага
            var h = if(h_ch == 0.0) (e/delta).pow(1.0/3.0) else h_ch

            var f1 = kutta2Const(syst, startVals, end, e, h)
            var f2 = kutta2Const(syst, startVals, end, e, h/2)

            val writer = FileWriter(File("src/main/resources/kutta2fullR.txt"))

            val rFull = MutableList<Double>(f1[0].size) {0.0}
            for(j in 0..<f1.size-1) {
                for (i in 0..<rFull.size) {
                    rFull[i] = (f2[2*j][i] - f1[j][i]) / (2.0.pow(2) - 1)
                }
                writer.write("${norm2(rFull.toDoubleArray())}\n")
            }

            while(norm2(rFull.toDoubleArray()) > e) {
                f1 = f2
                h /= 2
                f2 = kutta2Const(syst, startVals, end, e, h/2)


                writer.flush()
                for(j in 0..<f1.size-1) {
                    for (i in 0..<rFull.size) {
                        rFull[i] = (f2[2*j][i] - f1[j][i]) / (2.0.pow(2) - 1)
                    }
                    writer.write("${norm2(rFull.toDoubleArray())}\n")
                }
            }

            writer.close()

            return f2
        }
        private fun kutta2Const(syst: (DoubleArray) -> DoubleArray, startVals: DoubleArray, end:Double, e: Double, h_ch: Double = 0.0): MutableList<DoubleArray> {
            val c = 0.1
            val b = mutableListOf(-4.0 ,5.0)
            val a = 0.1

            val delta = (1.0/max(startVals[0], end)).pow(3) + norm2(syst(startVals)).pow(3) //выбор оптимального шага
            val h = if(h_ch == 0.0) (e/delta).pow(1.0/3.0) else h_ch

            val alls = mutableListOf(startVals)

            val nextVals = MutableList<Double>(startVals.size) {0.0}
            var cPos = startVals[0]

            while(cPos < end) {
                cPos += h
                if(cPos > end) cPos = end
                nextVals[0] = cPos //первым идет х - координата


                val k1 = multArr(syst(alls.last()), h)

                val coefVals = alls.last().toMutableList() //подготовим аргументы для вычисления k2
                coefVals[0] += c*h
                for(i in 1..< coefVals.size) {
                    coefVals[i] += a*k1[i-1]
                }
                val k2 = multArr(syst(coefVals.toDoubleArray()), h) //вычислили вектор второго коэфициента

                for(i in 1 ..< nextVals.size) {
                    nextVals[i] = alls.last()[i] + b[0]*k1[i-1] + b[1]*k2[i-1]
                }

                alls.add(nextVals.toDoubleArray())

            }

            return alls

        }

        fun kutta3ConstCheck(syst: (DoubleArray) -> DoubleArray, startVals: DoubleArray, end:Double, e: Double, h_ch: Double = 0.0): MutableList<DoubleArray> {

            val delta = (1.0/max(startVals[0], end)).pow(3) + norm2(syst(startVals)).pow(3) //выбор оптимального шага
            var h = if(h_ch == 0.0) (e/delta).pow(1.0/3.0) else h_ch

            var f1 = kutta3Const(syst, startVals, end, e, h)
            var f2 = kutta3Const(syst, startVals, end, e, h/2)

            val writer = FileWriter(File("src/main/resources/kutta3fullR.txt"))

            val rFull = MutableList<Double>(f1[0].size) {0.0}
            for(j in 0..<f1.size-1) {
                for (i in 0..<rFull.size) {
                    rFull[i] = (f2[2*j][i] - f1[j][i]) / (2.0.pow(2) - 1)
                }
                writer.write("${norm2(rFull.toDoubleArray())}\n")
            }

            while(norm2(rFull.toDoubleArray()) > e) {
                f1 = f2
                h /= 2
                f2 = kutta2Const(syst, startVals, end, e, h/2)

                writer.flush()
                for(j in 0..<f1.size-1) {
                    for (i in 0..<rFull.size) {
                        rFull[i] = (f2[2*j][i] - f1[j][i]) / (2.0.pow(2) - 1)
                    }
                    writer.write("${norm2(rFull.toDoubleArray())}\n")
                }
            }

            writer.close()

            return f2
        }
        private fun kutta3Const(syst: (DoubleArray) -> DoubleArray, startVals: DoubleArray, end:Double, e: Double, h_ch: Double = 0.0): MutableList<DoubleArray> {

            val delta = (1.0/max(startVals[0], end)).pow(3) + norm2(syst(startVals)).pow(3) //выбор оптимального шага
            val h = if(h_ch == 0.0) (e/delta).pow(1.0/3.0) else h_ch

            val alls = mutableListOf(startVals)

            val nextVals = MutableList<Double>(startVals.size) {0.0}
            var cPos = startVals[0]

            while(cPos < end) {
                cPos += h
                if(cPos > end) cPos = end
                nextVals[0] = cPos //первым идет х - координата


                val k1 = multArr(syst(alls.last()), h)

                var coefVals = alls.last().toMutableList() //подготовим аргументы для вычисления k2
                coefVals[0] += 0.5*h
                for(i in 1..< coefVals.size) {
                    coefVals[i] += 0.5*k1[i-1]
                }
                val k2 = multArr(syst(coefVals.toDoubleArray()), h) //вычислили вектор второго коэфициента

                coefVals = alls.last().toMutableList()
                coefVals[0] += h
                for(i in 1..< coefVals.size) {
                    coefVals[i] += -1*k1[i-1] + 2*k2[i-1]
                }
                val k3 = multArr(syst(coefVals.toDoubleArray()), h)//вычислили вектор третьего коэфициента

                for(i in 1 ..< nextVals.size) {
                    nextVals[i] = alls.last()[i] + 1.0/6.0*(k1[i-1] + 4*k2[i-1] + k3[i-1])
                }

                alls.add(nextVals.toDoubleArray())

            }

            return alls
        }


        fun kutta3Dynamic(syst: (DoubleArray) -> DoubleArray, startVals: DoubleArray, end:Double, e: Double, h_ch: Double = 0.0): MutableList<DoubleArray> {

            val delta = (1.0/max(startVals[0], end)).pow(3) + norm2(syst(startVals)).pow(3) //выбор оптимального шага
            var h = if(h_ch == 0.0) (e/delta).pow(1.0/3.0) else h_ch

            val alls = mutableListOf(startVals)
            val disp = mutableListOf(0.0)


            val nextVals = MutableList<Double>(startVals.size) {0.0}
            var cX = startVals[0]

            while(cX < end) {

                cX += h
                if(cX > end) {
                    h = cX - end
                    cX = end
                }
                nextVals[0] = cX //первым идет х - координата

                //////////////////////////////////////////////////////////// тут для выбранной точки считается приближение (с шагом h)
                var k1 = multArr(syst(alls.last()), h)

                var coefVals = alls.last().toMutableList() //подготовим аргументы для вычисления k2
                coefVals[0] += 0.5*h
                for(i in 1..< coefVals.size) {
                    coefVals[i] += 0.5*k1[i-1]
                }
                var k2 = multArr(syst(coefVals.toDoubleArray()), h) //вычислили вектор второго коэфициента

                coefVals = alls.last().toMutableList()
                coefVals[0] += h
                for(i in 1..< coefVals.size) {
                    coefVals[i] += -1*k1[i-1] + 2*k2[i-1]
                }
                var k3 = multArr(syst(coefVals.toDoubleArray()), h)

                for(i in 1 ..< nextVals.size) { //задаем новый набор точек - приближенного решения задачи Коши с шагом h
                    nextVals[i] = alls.last()[i] + 1.0/6.0*(k1[i-1] + 4*k2[i-1] + k3[i-1])
                }

                ///////////////////////////////////////////////////////// тут вычисляем значения приближений для половинного шага, определяем локальную погрешность
                val buffVals = MutableList<Double>(startVals.size) {0.0}
                buffVals[0] = cX - h/2

                k1 = multArr(syst(alls.last()), h/2)

                coefVals = alls.last().toMutableList() //подготовим аргументы для вычисления k2
                coefVals[0] += 0.5*h/2
                for(i in 1..< coefVals.size) {
                    coefVals[i] += 0.5*k1[i-1]
                }
                k2 = multArr(syst(coefVals.toDoubleArray()), h/2)

                coefVals = alls.last().toMutableList()
                coefVals[0] += h
                for(i in 1..< coefVals.size) {
                    coefVals[i] += -1*k1[i-1] + 2*k2[i-1]
                }
                k3 = multArr(syst(coefVals.toDoubleArray()), h)

                for(i in 1 ..< buffVals.size) { //задаем новый набор точек - приближенного решения задачи Коши с шагом h/2
                    buffVals[i] = alls.last()[i] + 1.0/6.0*(k1[i-1] + 4*k2[i-1] + k3[i-1])
                }
                //-------------------------------------------------------------// Вычислили в точке x0+ h/2. Теперь вычислим вторую точку (В x0+ h) c половинным шагом
                k1 = multArr(syst(buffVals.toDoubleArray()), h/2)

                coefVals = buffVals //подготовим аргументы для вычисления k2
                coefVals[0] += 0.5*h/2
                for(i in 1..< coefVals.size) {
                    coefVals[i] += 0.5*k1[i-1]
                }
                k2 = multArr(syst(coefVals.toDoubleArray()), h/2)

                coefVals = buffVals //подготовим аргументы для вычисления k2
                coefVals[0] += h/2
                for(i in 1..< coefVals.size) {
                    coefVals[i] +=  -1*k1[i-1] + 2*k2[i-1]
                }
                k3 = multArr(syst(coefVals.toDoubleArray()), h/2)


                for(i in 1 ..< buffVals.size) { //задаем новый набор точек - приближенного решения задачи Коши с шагом h/2
                    buffVals[i] = alls.last()[i] + 1.0/6.0*(k1[i-1] + 4*k2[i-1] + k3[i-1])
                }
                buffVals[0] = cX
                //-----------------------------------------------------------//в buffVals хранится занчение точки для нахождения локальной погрешности
                val cDisp = norm2(subtractVects(buffVals.toDoubleArray(), nextVals.toDoubleArray()))/(1.0 - 2.0.pow(-3))
                /////////////////////////////////////////////////////////


                if(cDisp < e) {
                    disp.add(cDisp)
                    alls.add(nextVals.toDoubleArray())
                    if(cDisp < e/2.0.pow(4)) {
                        h *= 2
                    }
                } else if(cDisp < e*2.0.pow(4)){
                    disp.add( norm2(subtractVects(buffVals.toDoubleArray(), nextVals.toDoubleArray()))/(2.0.pow(3) - 1.0))
                    alls.add(buffVals.toDoubleArray())
                    h /= 2
                } else {
                    cX -= h
                    h /= 2
                }

            }

            val w = FileWriter(File("src/main/resources/kutta3locDisp"))
            for(i in disp) {
                w.write("$i ")
            }
            w.close()

            return alls
        }
        fun kutta2Dynamic(syst: (DoubleArray) -> DoubleArray, startVals: DoubleArray, end:Double, e: Double, h_ch: Double = 0.0): MutableList<DoubleArray> {
            val c = 0.1
            val b = mutableListOf(-4.0 ,5.0)
            val a = 0.1

            val delta = (1.0/max(startVals[0], end)).pow(3) + norm2(syst(startVals)).pow(3) //выбор оптимального шага
            var h = if(h_ch == 0.0) (e/delta).pow(1.0/3.0) else h_ch

            val alls = mutableListOf(startVals)
            val disp = mutableListOf(0.0)


            val nextVals = MutableList<Double>(startVals.size) {0.0}
            var cX = startVals[0]

            while(cX < end) {

                cX += h
                if(cX > end) {
                    h = cX - end
                    cX = end
                }
                nextVals[0] = cX //первым идет х - координата

                //////////////////////////////////////////////////////////// тут для выбранной точки считается приближение (с шагом h)
                var k1 = multArr(syst(alls.last()), h)

                var coefVals = alls.last().toMutableList() //подготовим аргументы для вычисления k2
                coefVals[0] += c*h
                for(i in 1..< coefVals.size) {
                    coefVals[i] += a*k1[i-1]
                }
                var k2 = multArr(syst(coefVals.toDoubleArray()), h) //вычислили вектор второго коэфициента

                for(i in 1 ..< nextVals.size) { //задаем новый набор точек - приближенного решения задачи Коши с шагом h
                    nextVals[i] = alls.last()[i] + b[0]*k1[i-1] + b[1]*k2[i-1]
                }

                ///////////////////////////////////////////////////////// тут вычисляем значения приближений для половинного шага, определяем локальную погрешность
                val buffVals = MutableList<Double>(startVals.size) {0.0}
                buffVals[0] = cX - h/2

                k1 = multArr(syst(alls.last()), h/2)

                coefVals = alls.last().toMutableList() //подготовим аргументы для вычисления k2
                coefVals[0] += c*h/2
                for(i in 1..< coefVals.size) {
                    coefVals[i] += a*k1[i-1]
                }
                k2 = multArr(syst(coefVals.toDoubleArray()), h/2)

                for(i in 1 ..< buffVals.size) { //задаем новый набор точек - приближенного решения задачи Коши с шагом h/2
                    buffVals[i] = alls.last()[i] + b[0]*k1[i-1] + b[1]*k2[i-1]
                }
                //-------------------------------------------------------------// Вычислили в точке x0+ h/2. Теперь вычислим вторую точку (В x0+ h) c половинным шагом
                k1 = multArr(syst(buffVals.toDoubleArray()), h/2)

                coefVals = buffVals //подготовим аргументы для вычисления k2
                coefVals[0] += c*h/2
                for(i in 1..< coefVals.size) {
                    coefVals[i] += a*k1[i-1]
                }
                k2 = multArr(syst(coefVals.toDoubleArray()), h/2)

                for(i in 1 ..< buffVals.size) { //задаем новый набор точек - приближенного решения задачи Коши с шагом h/2
                    buffVals[i] = alls.last()[i] + b[0]*k1[i-1] + b[1]*k2[i-1]
                }
                buffVals[0] = cX
                //-----------------------------------------------------------//в buffVals хранится занчение точки для нахождения локальной погрешности
                val cDisp = norm2(subtractVects(buffVals.toDoubleArray(), nextVals.toDoubleArray()))/(1.0 - 2.0.pow(-2))
                /////////////////////////////////////////////////////////


                if(cDisp < e) {
                    disp.add(cDisp)
                    alls.add(nextVals.toDoubleArray())
                    if(cDisp < e/2.0.pow(3)) {
                        h *= 2
                    }
                } else if(cDisp < e*2.0.pow(3)){
                    disp.add( norm2(subtractVects(buffVals.toDoubleArray(), nextVals.toDoubleArray()))/(2.0.pow(2) - 1.0))
                    alls.add(buffVals.toDoubleArray())
                    h /= 2
                } else {
                    cX -= h
                    h /= 2
                }

            }
            val w = FileWriter(File("src/main/resources/kutta2locDisp"))
            for(i in disp) {
                w.write("$i ")
            }
            w.close()

            return alls
        }

        private fun multArr(vect: DoubleArray, h: Double): DoubleArray {
            val resArr = vect.toMutableList()
            for (i in vect.indices) {
                resArr[i] *= h
            }
            return resArr.toDoubleArray()
        }
        private fun subtractVects(vect1: DoubleArray, vect2: DoubleArray): DoubleArray {
            val res = DoubleArray(vect1.size) {0.0}
            for(i in 0..< vect1.size) {
                res[i] = vect1[i] - vect2[i]
            }
            return res
        }
        private fun norm2( vect: DoubleArray): Double {
            var res = 0.0
            for(i in vect) {
                res += i*i
            }
            res = sqrt(res)

            return res
        }
    }
}