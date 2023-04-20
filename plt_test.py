import matplotlib.pyplot as plt
save_path = '/Users/josephko/research/ice_renders/20230413'
plt.plot([0, 1, 2, 3, 4], [0, 3, 5, 9, 11])
file_path = save_path + '/' + 'test.png'
plt.savefig(file_path)
plt.show()