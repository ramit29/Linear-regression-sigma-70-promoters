import web
from Bio import motifs
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
#from bs4 import BeautifulSoup
import time




def make_text(string):
    return string

urls = ('/', 'tutorial')
render = web.template.render('templates/')

app = web.application(urls, globals())

my_form = web.form.Form(
                web.form.Textbox('', class_='textfield', id='textfield'),
                )

class tutorial:
    def GET(self):
        form = my_form()
        return render.tutorial(form, "")
        
    def POST(self):
        form = my_form()
        form.validates()
        s = form.value['textfield']
        print s
        
        s1 = s[0:6]
        s2 = s[7:13]
        print s1
        print s2
        #j = Seq(s);
        #print j
       
        # -10 Samples

        instances = [Seq("TACAGT"),
        Seq("TATTAT"),
        Seq("TACTGT"),
        Seq("GATTAT"),
        Seq("TATTGT"),
        Seq("TACTAT"),
        Seq("TATAGT"),
        Seq("TATTAT"),
        Seq("TATAAT"),
        Seq("GACTGT"),
        Seq("TACAAT"),
        Seq("TATAGT"),
        Seq("GATTAT"),
        Seq("GATTAT"),
        Seq("TACAAT"),
        Seq("TACAAT"),
        Seq("GACTAT"),
        Seq("GATTGT"),
        Seq("TATTGT")
        ]

        """
        instances1 = [Seq("TTGACG"),
        Seq("TTTACA"),
        Seq("TTGACA"),
        Seq("CTGATA"),
        Seq("TTGACA"),
        Seq("TTTACG"),
        Seq(s)
        ]
        print instances1
        """
        instancesP = [Seq("TACAGT"),
        Seq("TATTAT"),
        Seq("TACTGT"),
        Seq("GATTAT"),
        Seq("TATTGT"),
        Seq("TACTAT"),
        Seq("TATAGT"),
        Seq("TATTAT"),
        Seq("TATAAT"),
        Seq("GACTGT"),
        Seq("TACAAT"),
        Seq("TATAGT"),
        Seq("GATTAT"),
        Seq("GATTAT"),
        Seq("TACAAT"),
        Seq("TACAAT"),
        Seq("GACTAT"),
        Seq("GATTGT"),
        Seq("TATTGT"),
        Seq(s1)
        ]

        m = motifs.create(instances)
        #print(m.counts);
        pwm = m.counts.normalize(pseudocounts=0.5)
        #print(pwm)
        pssm = pwm.log_odds()
        #print(pssm)
        #return pssm
        result = [[0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0]
        ] 

        mP = motifs.create(instancesP)
        print(mP.counts)
        pwmP = mP.counts.normalize(pseudocounts=0.5)
        print(pwmP)
        pssmP = pwmP.log_odds()
        print(pssmP)
        p,o,l,k,m,n = str(s1)
        print "bro"
        print p
        print o
        print l
        print k
        print m
        print n
        #return pssm
        resultP = pssmP[p,0] + pssmP[o,1] + pssmP[l,2] + pssmP[k,3] + pssmP[m,4] + pssmP[n,5]
        print resultP


        def calculateX(a,b,c,d,e,f,x):
            result[x][0] = pssm[a,0] + pssm[b,1] + pssm[c,2] + pssm[d,3] + pssm[e,4] + pssm[f,5]
            return result[x][0] 
            calculateX('T', 'A', 'C', 'A', 'G', 'T',0)    
            calculateX('T', 'A', 'T', 'T', 'A', 'T',1)
            calculateX('T', 'A', 'C', 'T', 'G', 'T',2)    
            calculateX('G', 'A', 'T', 'T', 'A', 'T',3)
            calculateX('T', 'A', 'T', 'T', 'G', 'T',4)    
            calculateX('T', 'A', 'C', 'T', 'A', 'T',5)
            calculateX('T', 'A', 'T', 'A', 'G', 'T',6)    
            calculateX('T', 'A', 'T', 'T', 'A', 'T',7)
            calculateX('T', 'A', 'T', 'A', 'A', 'T',8)    
            calculateX('G', 'A', 'C', 'T', 'G', 'T',9)
            calculateX('T', 'A', 'C', 'A', 'A', 'T',10)
            calculateX('T', 'A', 'T', 'A', 'G', 'T',11)
            calculateX('G', 'A', 'T', 'T', 'A', 'T',12)    
            calculateX('G', 'A', 'T', 'T', 'A', 'T',13)
            calculateX('T', 'A', 'C', 'A', 'A', 'T',14)    
            calculateX('T', 'A', 'C', 'A', 'A', 'T',15)
            calculateX('G', 'A', 'C', 'T', 'A', 'T',16)    
            calculateX('G', 'A', 'T', 'T', 'G', 'T',17)
            calculateX('T', 'A', 'T', 'T', 'G', 'T',18) 
        print result  
        outputResult = [[1],
        [0.7],
        [0.86],
        [0.01],
        [0.72],
        [0.24],
        [0.47],
        [0.36],
        [0.51],
        [0.04],
        [0.33],
        [0.58],
        [0],
        [0.01],
        [0.1],
        [0.15],
        [0.16],
        [0.06],
        [0.56]
        ]
        #print result
        #print outputResult        
        #return outputResult
        #return make_text(s)


        # -35 Samples

        instances2 = [Seq("TTGACG"),
        Seq("TTTACA"),
        Seq("TTGACA"),
        Seq("CTGATA"),
        Seq("TTGACA"),
        Seq("TTTACG"),
        Seq("TTTACG"),
        Seq("TTTACG"),
        Seq("CTGACA"),
        Seq("TTTACA"),
        Seq("TTTACG"),
        Seq("TTGACG"),
        Seq("CTGATA"),
        Seq("CTGATG"),
        Seq("TTTATG"),
        Seq("TTTATA"),
        Seq("TTGACA"),
        Seq("TTGACA"),
        Seq("TTGACG")
        ]

        instancesP2 = [Seq("TTGACG"),
        Seq("TTTACA"),
        Seq("TTGACA"),
        Seq("CTGATA"),
        Seq("TTGACA"),
        Seq("TTTACG"),
        Seq("TTTACG"),
        Seq("TTTACG"),
        Seq("CTGACA"),
        Seq("TTTACA"),
        Seq("TTTACG"),
        Seq("TTGACG"),
        Seq("CTGATA"),
        Seq("CTGATG"),
        Seq("TTTATG"),
        Seq("TTTATA"),
        Seq("TTGACA"),
        Seq("TTGACA"),
        Seq("TTGACG"),
        Seq(s2)
        ]
        m2 = motifs.create(instances2)
        print(m2.counts)
        pwm2 = m2.counts.normalize(pseudocounts=0.5)
        print(pwm2)
        pssm2 = pwm2.log_odds()
        print(pssm2)
        result2 = [[0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0],
        [0]
        ] 
        
        mP2 = motifs.create(instancesP2)
        print(mP2.counts)
        pwmP2 = mP2.counts.normalize(pseudocounts=0.5)
        print(pwmP2)
        pssmP2 = pwmP2.log_odds()
        print(pssmP2)
        p2,o2,l2,k2,m2,n2 = str(s2)
        print "bro"
        print p2
        print o2
        print l2
        print k2
        print m2
        print n2
        #return pssm2
        resultP2 = pssmP[p2,0] + pssmP[o2,1] + pssmP[l2,2] + pssmP[k2,3] + pssmP[m2,4] + pssmP[n2,5]
        print resultP2


        def calculateX2(a,b,c,d,e,f,x):
            result2[x][0] = pssm[a,0] + pssm[b,1] + pssm[c,2] + pssm[d,3] + pssm[e,4] + pssm[f,5]
            return result[x][0] 
        calculateX2('T', 'T', 'G', 'A', 'C', 'G',0)    
        calculateX2('T', 'T', 'T', 'A', 'C', 'A',1)
        calculateX2('T', 'T', 'G', 'A', 'C', 'A',2)    
        calculateX2('C', 'T', 'G', 'A', 'T', 'A',3)
        calculateX2('T', 'T', 'G', 'A', 'C', 'A',4)    
        calculateX2('T', 'T', 'T', 'A', 'C', 'G',5)
        calculateX2('T', 'T', 'T', 'A', 'C', 'G',6)    
        calculateX2('T', 'T', 'T', 'A', 'C', 'G',7)
        calculateX2('C', 'T', 'G', 'A', 'C', 'A',8)    
        calculateX2('T', 'T', 'T', 'A', 'C', 'A',9)
        calculateX2('T', 'T', 'T', 'A', 'C', 'G',10)
        calculateX2('T', 'T', 'G', 'A', 'C', 'G',11)
        calculateX2('C', 'T', 'G', 'A', 'T', 'A',12)    
        calculateX2('C', 'T', 'G', 'A', 'T', 'A',13)
        calculateX2('T', 'T', 'T', 'A', 'T', 'G',14)    
        calculateX2('T', 'T', 'T', 'A', 'T', 'A',15)
        calculateX2('T', 'T', 'G', 'A', 'C', 'A',16)    
        calculateX2('T', 'T', 'G', 'A', 'C', 'A',17)
        calculateX2('T', 'T', 'G', 'A', 'C', 'G',18)    
        print result2

        outputResult2 = [[1],
        [0.7],
        [0.86],
        [0.01],
        [0.72],
        [0.24],
        [0.47],
        [0.36],
        [0.51],
        [0.04],
        [0.33],
        [0.58],
        [0],
        [0.01],
        [0.1],
        [0.15],
        [0.16],
        [0.06],
        [0.56]
        ]
        #print outputResult2 

        #Combined output result

        theoutputResult = [[1],
        [0.7],
        [0.86],
        [0.01],
        [0.72],
        [0.24],
        [0.47],
        [0.36],
        [0.51],
        [0.04],
        [0.33],
        [0.58],
        [0],
        [0.01],
        [0.1],
        [0.15],
        [0.16],
        [0.06],
        [0.56],
        [1],
        [0.7],
        [0.86],
        [0.01],
        [0.72],
        [0.24],
        [0.47],
        [0.36],
        [0.51],
        [0.04],
        [0.33],
        [0.58],
        [0],
        [0.01],
        [0.1],
        [0.15],
        [0.16],
        [0.06],
        [0.56]
        ]

        inputResult = result + result2
        print "bro"
        print len(inputResult)
        print len(theoutputResult)

        model = LinearRegression() 
        model.fit(inputResult,theoutputResult)
        x2 = np.linspace(5.5,9.5)
        y2 = model.predict(x2[:,None])
        #plt.scatter(result,outputResult,color='#26a69a')
        plt.scatter(result,outputPartialResult,color='blue')
        plt.scatter(result2,outputPartialResult,color='green')
        plt.plot(x2,y2,'r')
        plt.legend(['Predicted','Observed'])
        plt.xlabel('Input', fontsize=16)
        plt.ylabel('Output', fontsize=16)
        #ax = plt.gca()
        #leg = ax.get_legend()
        #leg.legendHandles[0].set_color('red')
        #leg.legendHandles[1].set_color('green')
        #Pridict y value for an x value
        prediction = model.predict(resultP)
        pre = str(prediction)
        fig1 = plt.gcf()
        #plt.show()
        fig1.savefig('books_read3.png', dpi=100)
        #source_code = """<span class="UserName"><a href="#">Martin Elias</a></span>"""
        #source_code = """<img src="books_read1.png">"""
        #soup = BeautifulSoup(source_code)
        #print soup.a.string
        data_uri = open('books_read3.png', 'rb').read().encode('base64').replace('\n', '')
        

        d1 = '<div class="row"><div class="col s5 offset-s1">'
        graph = '<h2>Regression Graph :</h2>'
        img_tag = '<img src="data:image/png;base64,%s">' % data_uri
        d2 = '</div>'
        d3 = '<div class="col s4 offset-s2">'
        start = '<h2>Key : </h2><br>'
        key1 = '<p>Observed : <i id="red" class="fa fa-long-arrow-right" aria-hidden="true"></i></p><br>'
        key2 = '<p>Pridicted : <i id = green class="fa fa-circle" aria-hidden="true"></i></p><br>'
        img_tag1 = '<p>Input :<span> %s </span></p><br>' % s 
        img_tag2 = '<p>X value of Input :<span> %s </span></p><br>' % resultP
        img_tag3 = '<p>Predicted Value :<span> %s <span/></p>' % pre
        return d1 + graph + img_tag + d2 + d3 + start + key1 + key2 + img_tag1 + img_tag2 + img_tag3 + d2 + d2
        #return plt.show()

   

if __name__ == '__main__':
    app.run()

