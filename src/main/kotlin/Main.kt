import java.io.File
import java.io.FileWriter
import kotlin.math.pow

//A = 1/10
//B = 1/20
//ksi = 1/10

var calls = 0
fun main(args: Array<String>) {


    val w2 = FileWriter(File("src/main/resources/kutta2.txt"))
    val w3 = FileWriter(File("src/main/resources/kutta3.txt"))

    val firstSpot = doubleArrayOf(0.0, 3.14/20.0, 3.14/10.0)
    val end = 3.14
    val e = 10.0.pow(-5.0)

    fun syst(y: DoubleArray): DoubleArray {
        calls ++
        //*помни, что первым элементом в массиве передается позиция х

        var result = mutableListOf(y[2]/10.0, -y[1]/20.0)
        return result.toDoubleArray()

    }
    //////////////////////////////////////////////////////////////////////////// до этого момента были вводные данные, сейчас же само выполнение интегрирования
    var eqData2 = koshiSolver.kutta2Dynamic(::syst, firstSpot, end, e)
    println("При выполнении kutta2 правая часть была вызвана $calls раз")
    calls = 0
    for(i in 0..<eqData2[0].size) {
        for(j in 0..<eqData2.size) {
            w2.write("${eqData2[j][i]} ")
        }
        w2.write("\n")
    }


    var eqData3 = koshiSolver.kutta3Dynamic(::syst, firstSpot, end, e)
    println("При выполнении kutta3 правая часть была вызвана $calls раз")
    calls = 0
    for(i in 0..<eqData3[0].size) {
        for(j in 0..<eqData3.size) {
            w3.write("${eqData3[j][i]} ")
        }
        w3.write("\n")
    }
//
//    calls = 0
//    for(i in 0..5) {
//        koshiSolver.kutta2Dynamic(::syst, firstSpot, end, 10.0.pow(-1*i))
//        print("${calls} ")
//        calls = 0
//    }
//    println()
//
//    calls = 0
//    for(i in 0..5) {
//        koshiSolver.kutta3Dynamic(::syst, firstSpot, end, 10.0.pow(-1*i))
//        print("${calls} ")
//        calls = 0
//    }
//

    w2.close()
    w3.close()
}