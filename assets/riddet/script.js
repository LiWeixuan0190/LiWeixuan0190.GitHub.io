var searchButton = document.getElementById("navbar-search-button");
searchButton.addEventListener("click", search);

var searchField = document.getElementById("navbar-search-input");
searchField.addEventListener("input", search);


function search()
{
		var searchValue = searchField.value;
		searchValue = searchValue.toLowerCase();
		searchValue = searchValue.trim();

		Array.from(document.getElementsByClassName("post")).forEach(function(item){
			if (item.textContent.toLowerCase().indexOf(searchValue) == -1)
			{
				item.style.display = "none";
			}
			if (item.textContent.toLowerCase().indexOf(searchValue) != -1)
			{
				item.style.display = "";
			}
		});
}

function createTwit()
{
	if(inputText.value == "" || inputAuthor.value == "")
	{
		alert("Please fill out both fields!");
		return;
	}

	var newTwit = document.createElement("article");
	newTwit.className = "twit";

	// Icon
	var icon_div = document.createElement("div");
	icon_div.className = "twit-icon";
	var icon = document.createElement("i");
	icon.className = "fa fa-bullhorn";
	icon_div.appendChild(icon);
	newTwit.appendChild(icon_div);

	// Text
	var textNode = document.createElement("p");
	textNode.className = "twit-text";
	var twitText = document.createTextNode(inputText.value);
	textNode.appendChild(twitText);

	// Author
	var authorNode = document.createElement("p");
	authorNode.className = "twit-attribution";
	var link = document.createElement("a");
	var authorText = document.createTextNode(inputAuthor.value);
	link.appendChild(authorText);
	link.href = "#";
	authorNode.appendChild(link);


	var content = document.createElement("div");
	content.className = "twit-content";
	content.appendChild(textNode);
	content.appendChild(authorNode);

	newTwit.appendChild(content);

	document.getElementsByClassName("twit-container")[0].appendChild(newTwit);

	hideModal();
}
