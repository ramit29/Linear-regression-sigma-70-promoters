import matplotlib.pyplot as plt

plt.plot([0, 1, 2, 3, 4], [0, 3, 5, 9, 11])

plt.xlabel('Months')
plt.ylabel('Books Read')
#plt.savefig()

fig1 = plt.gcf()
plt.show()
fig1.savefig('books_read1.png', dpi=100)	